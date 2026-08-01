#pragma once
#include "cryptoTools/Common/Defines.h"
#include "cryptoTools/Common/block.h"
#include "cryptoTools/Common/BitVector.h"
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Common/Timer.h"
#include "cryptoTools/Network/Channel.h"
#include "coproto/coproto.h"

#include "seal/seal.h"

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
        // Number of layers actually transmitted. It is a public function of
        // (n, heNumSlots, m, span_blocks) and must not depend on the sender's
        // set: the realized layer count does, so sending it would leak. Left
        // at 0 the parties derive the certified budget of
        // certifiedLayerBudget() at init; set it explicitly only to explore
        // the (uncertified) empirical regime.
        u32 layerBudget = 0;
        // Statistical security parameter for the layer budget.
        u32 layerBudgetLambda = 40;
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

    // Combinatorial core of the layer sequencing: partition items (bin, block)
    // into layers such that each layer holds at most one item per bin and the
    // occupied blocks of a layer span at most spanBlocks consecutive blocks.
    // Processes blocks left to right, assigns each item to the compatible
    // layer with the smallest anchor block (earliest-expiring first), and
    // opens a new layer anchored at the current block when none fits. This is
    // exactly the optimal algorithm for the underlying interval-covering LP;
    // sequencingLowerBound() certifies per-instance optimality.
    // Returns the layer count and fills itemToLayer / layerMinBlock /
    // layerMaxBlock (resized internally).
    uint32_t sequenceLayers(
        const std::vector<uint32_t> &itemBin,
        const std::vector<uint32_t> &itemBlock,
        uint32_t numSlots, uint32_t spanBlocks,
        std::vector<uint32_t> &itemToLayer,
        std::vector<uint32_t> &layerMinBlock,
        std::vector<uint32_t> &layerMaxBlock);

    // Exact lower bound on the achievable layer count, via the LP dual of the
    // sequencing problem: the maximum over families of block intervals
    // J_1..J_K whose spanBlocks-extensions are pairwise disjoint of
    // sum_k (max per-bin item count inside J_k). Every valid layer partition
    // needs at least this many layers (each layer's anchor block lies in at
    // most one extension), so sequenceLayers(...) == sequencingLowerBound(...)
    // certifies that both are optimal for the given instance.
    uint64_t sequencingLowerBound(
        const std::vector<uint32_t> &itemBin,
        const std::vector<uint32_t> &itemBlock,
        uint32_t numSlots, uint32_t spanBlocks);

    // Public layer budget: a number of layers that the optimal sequencing of
    // n uniformly placed items fits into, except with probability at most
    // 2^-lambda. The protocol always transmits exactly this many layers (real
    // ones plus padding), so the realized layer count -- a function of the
    // sender's set -- never reaches the receiver.
    //
    // Derivation. A layer is a window of spanBlocks consecutive blocks with
    // capacity one per bin, so a multiset of window starts admits all items
    // iff Hall's condition holds for every block interval J of length g:
    // the windows overlapping J must number at least the largest per-bin item
    // count inside J. Laying windows down at a uniform density rho satisfies
    // every constraint once rho*(g+spanBlocks) - 1 >= D(g) for all g, where
    // D(g) bounds that per-bin count. D(g) is the exact upper binomial
    // quantile of Bin(n, g/positionRange) at the union-adjusted level
    // 2^-lambda / (numSlots * b^2) -- one budget share per (bin, position,
    // length) constraint -- and the resulting window count is
    // ceil(rho * (b + spanBlocks)).
    uint32_t certifiedLayerBudget(
        uint64_t n, uint32_t numSlots, uint32_t numBlocks,
        uint64_t positionRange, uint32_t spanBlocks, uint32_t lambda = 40);

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
        std::vector<std::vector<uint32_t>> mLayerBins;
        std::vector<uint32_t> mLayerMinBlock;
        std::vector<uint32_t> mLayerMaxBlock;
        std::vector<uint8_t> mLayerIsPadding;

        // Slot -> item index of the item sitting there, UINT32_MAX if empty,
        // in layer-major order. Lets the OT-based extensions (PSU, card-sum)
        // address the sender's payloads by slot.
        std::vector<uint32_t> mSlotToItem;
        // One mask per slot, layer-major: mask[lay * mNumSlots + bin]. The
        // equality runs over every slot of the L x H layout, so the occupancy
        // is never transmitted -- that disclosure is exactly the
        // layout-matching leak of the KKLS follow-up note (the receiver's band
        // hash is public, so it can predict its own items' slot columns).
        std::vector<uint64_t> mMasks;
        std::vector<seal::Plaintext> ptxts_mask;

        std::vector<std::vector<seal::Plaintext>> ptxts_diags;

        uint64_t mOTeBatchSize = 1ull << 19;

    public:
        void sequencing(const std::vector<uint32_t>& start_pos_spacing);

        void init(
            uint32_t n, uint32_t nReceiver,
            sspmtParams ssParams, oc::block seed = oc::OneBlock);

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

        uint64_t mOTeBatchSize = 1ull << 19;

    public:
        void init(
            uint32_t n, uint32_t nSender,
            sspmtParams ssParams, oc::block seed = oc::ZeroBlock);

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
