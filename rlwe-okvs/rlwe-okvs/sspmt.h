#pragma once
#include "cryptoTools/Common/Defines.h"
#include "cryptoTools/Common/block.h"
#include "cryptoTools/Common/BitVector.h"
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Common/Timer.h"
#include "cryptoTools/Network/Channel.h"
#include "coproto/coproto.h"

#include "seal/seal.h"

#include "../GMW/Gmw.h"

using namespace std;
using namespace seal;
using namespace oc;

namespace rlweOkvs
{
    using Proto = coproto::task<>;
    using Socket = coproto::Socket;

    struct sspmtParams {
        u32 heNumSlots = 1 << 13;
        std::vector<int> heCoeffModulus;
        u32 hePlainModulusBits;
        u32 bandWidth;
        u32 span_blocks;
        double bandExpansion;
        // Number of layers actually transmitted; see resolveLayerBudget().
        // Left at 0 the parties derive the certified budget at init; set it
        // explicitly only to explore the (uncertified) empirical regime.
        u32 layerBudget = 0;
        // Statistical security parameter for the layer budget.
        u32 layerBudgetLambda = 40;

        // The number of layers to transmit for a sender set of size n: the
        // explicit layerBudget when set, otherwise a count that the optimal
        // sequencing of n uniformly placed items fits into except with
        // probability at most 2^-layerBudgetLambda. It is a function of the
        // public parameters alone, so both parties derive the same value
        // without communicating -- the realized layer count depends on the
        // sender's set, so sending that would leak.
        //
        // Derivation. A layer is a window of span_blocks consecutive blocks
        // with capacity one per bin, so a multiset of window starts admits all
        // items iff Hall's condition holds for every block interval J of
        // length g: the windows overlapping J must number at least the largest
        // per-bin item count inside J. Laying windows down at a uniform
        // density rho satisfies every constraint once
        // rho*(g+span_blocks) - 1 >= D(g) for all g, where D(g) bounds that
        // per-bin count. D(g) is the exact upper binomial quantile of
        // Bin(n, g/positionRange) at the union-adjusted level
        // 2^-lambda / (heNumSlots * b^2) -- one budget share per (bin,
        // position, length) constraint -- and the window count is then
        // ceil(rho * (b + span_blocks)).
        u32 resolveLayerBudget(u64 n) const;
        void initialize(int n) {
            switch(n){
                case (1ull << 16):
                    heCoeffModulus = {50, 58, 60, 50};
                    hePlainModulusBits = 56;
                    bandWidth = 28;
                    bandExpansion = 2.3;
                    span_blocks = 9;
                    break;
                case (1ull << 18):
                    heCoeffModulus = {54, 58, 60, 50};
                    hePlainModulusBits = 58;
                    bandWidth = 28;
                    bandExpansion = 2.2;
                    span_blocks = 13;
                    break;
                case (1ull << 20):
                    heCoeffModulus = {58, 58, 60, 50};
                    hePlainModulusBits = 60;
                    bandWidth = 31;
                    bandExpansion = 2.1;
                    span_blocks = 20;
                    break;
                case (1ull << 22):
                    heCoeffModulus = {60, 60, 60, 50};
                    hePlainModulusBits = 60;
                    bandWidth = 45;
                    bandExpansion = 1.7;
                    span_blocks = 30;
                    break;
                default:
                    heCoeffModulus = {58, 58, 60, 50};
                    hePlainModulusBits = 60;
                    bandWidth = 53;
                    bandExpansion = 1.5;
                    span_blocks = 20;
                    break;
            }
        }
    };

    class SspmtSender: public oc::TimerAdapter
    {

        oc::PRNG mPrng;
        Modulus mModulus;
        uint64_t mNumSlots;
        unique_ptr<Evaluator> mEvaluator;
        unique_ptr<BatchEncoder> mBatchEncoder;
        shared_ptr<SEALContext> mContext;

        uint32_t mN, mNreceiver, mM, mW, mNumBatch, mWrap;
        uint32_t mNumLayers;      // == mLayerBudget; padded, never data dependent
        uint32_t mNumRealLayers;  // sequencing output; sender-private
        uint32_t mLayerBudget;
        uint32_t mSpanBlocks;
        std::vector<uint32_t> mItemToLayerIdx;
        std::vector<uint32_t> mItemToBlockIdx;
        std::vector<uint32_t> mLayerMinBlock;
        std::vector<uint32_t> mLayerMaxBlock;
        std::vector<uint8_t> mLayerIsPadding;

        // The layout itself: slot -> index of the item sitting there,
        // UINT32_MAX if empty, layer-major (slot = layer * mNumSlots + bin).
        // Also what the OT-based extensions (PSU, card-sum) use to address the
        // sender's payloads by slot.
        std::vector<uint32_t> mSlotToItem;
        // One mask per slot, layer-major: mask[lay * mNumSlots + bin]. The
        // equality runs over every slot of the L x H layout, so the occupancy
        // is never transmitted -- that disclosure is exactly the
        // layout-matching leak of the KKLS follow-up note (the receiver's band
        // hash is public, so it can predict its own items' slot columns).
        std::vector<uint64_t> mMasks;
        std::vector<seal::Plaintext> ptxts_mask;

        std::vector<std::vector<seal::Plaintext>> ptxts_diags;

        // The equality GMW. Its size is public (the layout is), so the
        // correlated randomness it needs can be produced before either set is
        // known -- see setup().
        volePSI::Gmw mGmw;
        bool mSetupDone = false;
        uint64_t mOTeBatchSize = 1ull << 19;

        void initGmw();

    public:
        // Combinatorial core of the sequencing: partition items (bin, block)
        // into layers such that each layer holds at most one item per bin and
        // the occupied blocks of a layer span at most spanBlocks consecutive
        // blocks. Processes blocks left to right, assigns each item to the
        // compatible layer with the smallest anchor block (earliest-expiring
        // first), and opens a new layer anchored at the current block when
        // none fits. It has matched the optimum on every instance measured
        // (Sequencing_test checks it against both the LP dual and
        // sequenceLayersOptimal), but that is an observation, not a proof --
        // sequencing() falls back to sequenceLayersOptimal if it ever
        // overshoots the budget. Returns the layer count and fills
        // itemToLayer / layerMinBlock / layerMaxBlock (resized internally).
        // Static so it can be exercised on its own, without HE state.
        static uint32_t sequenceLayers(
            const std::vector<uint32_t>& itemBin,
            const std::vector<uint32_t>& itemBlock,
            uint32_t numSlots, uint32_t spanBlocks,
            std::vector<uint32_t>& itemToLayer,
            std::vector<uint32_t>& layerMinBlock,
            std::vector<uint32_t>& layerMaxBlock);

        // Minimum-layer sequencing, same interface as sequenceLayers.
        //
        // A layer is a window of spanBlocks consecutive blocks with capacity
        // one per bin, so by Hall's condition (which for these convex
        // neighbourhoods collapses to block intervals) a multiset of window
        // starts admits all items iff every interval J gets at least
        // M(J) = max_bin (items of that bin inside J) windows overlapping it.
        // Those supports are runs of consecutive starts, so the covering
        // program is an interval program: sweeping the right endpoint and
        // filling each violated constraint at the latest position that can
        // serve it (the current block) is optimal by the standard exchange.
        // Items are then assigned per bin, leftmost admissible window first,
        // which is optimal for convex bipartite matching.
        //
        // Costs O(n*b + b^2) against the greedy's O(n*L), so it is the slow
        // path: sequencing() runs the greedy and only falls back here if that
        // overshoots the budget.
        static uint32_t sequenceLayersOptimal(
            const std::vector<uint32_t>& itemBin,
            const std::vector<uint32_t>& itemBlock,
            uint32_t numSlots, uint32_t numBlocks, uint32_t spanBlocks,
            std::vector<uint32_t>& itemToLayer,
            std::vector<uint32_t>& layerMinBlock,
            std::vector<uint32_t>& layerMaxBlock);

        // Runs sequenceLayers on the band start positions, then pads the
        // result out to the public layer budget and builds the layout.
        void sequencing(const std::vector<uint32_t>& start_pos_spacing);

        void init(
            uint32_t n, uint32_t nReceiver,
            sspmtParams ssParams, oc::block seed = oc::OneBlock);

        // Offline phase: generate the GMW triples. Depends only on the public
        // parameters, so it can run at any point before run(); run() falls
        // back to generating them inline if it was skipped.
        Proto setup(Socket& chl);

        u32 getNumLayers() {return mNumLayers;};
        // Number of equality instances run: the whole L x H layout.
        u64 getLayoutSize() {return (u64)mNumLayers * mNumSlots;};
        // Slot -> item index (UINT32_MAX when the slot holds no item).
        const std::vector<uint32_t>& getSlotToItem() const {return mSlotToItem;};

        Proto run(
            const std::vector<oc::block> &Y,
            oc::BitVector &results,
            Socket &chl);

        void preprocess(
            const std::vector<oc::block> &Y);

        Proto recv_encoded_chunks(
            std::vector<std::vector<seal::Ciphertext>> &encoded_in_he,
            Socket &chl);

        Proto send_decoded_chunks(
            const std::vector<std::vector<seal::Ciphertext>> &encoded_in_he,
            Socket &chl);

        void encrypted_decode(
            const std::vector<std::vector<seal::Ciphertext>> &encoded_in_he,
            std::vector<seal::Ciphertext> &decoded_in_he,
            uint32_t layerBegin,
            uint32_t layerEnd);
    };

    class SspmtReceiver: public oc::TimerAdapter
    {
        oc::PRNG mPrng;
        Modulus mModulus;
        uint64_t mNumSlots;
        shared_ptr<SEALContext> mContext;
        unique_ptr<Encryptor> mEncryptor;
        unique_ptr<BatchEncoder> mBatchEncoder;
        unique_ptr<Decryptor> mDecryptor;

        uint32_t mN, mNsender, mM, mW, mNumBatch, mWrap;
        uint32_t mLayerBudget;

        uint64_t mIndicatorStr;

        volePSI::Gmw mGmw;
        bool mSetupDone = false;
        uint64_t mOTeBatchSize = 1ull << 19;

        void initGmw();

    public:
        void init(
            uint32_t n, uint32_t nSender,
            sspmtParams ssParams, oc::block seed = oc::ZeroBlock);

        // Offline phase; see SspmtSender::setup().
        Proto setup(Socket& chl);

        u32 getNumLayers() {return mLayerBudget;};
        u64 getLayoutSize() {return (u64)mLayerBudget * mNumSlots;};

        Proto run(
            const std::vector<oc::block> &X,
            oc::BitVector &results,
            Socket &chl);

        Proto send_encoded_chunks(
            const std::vector<oc::block> &X,
            Socket &chl);

        Proto recv_decoded_chunks(
            std::vector<seal::Ciphertext> &decoded_in_he,
            Socket &chl);
    };
}
