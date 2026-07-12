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
        uint32_t mW_seq;
        uint32_t mNumLayers;
        uint32_t mSpanBlocks;
        std::vector<uint32_t> mItemToLayerIdx;
        std::vector<uint32_t> mItemToBlockIdx;
        std::vector<std::vector<uint32_t>> mLayerBins;
        std::vector<uint32_t> mLayerMinBlock;
        std::vector<uint32_t> mLayerMaxBlock;
        
        std::vector<uint32_t> last_layer_per_bin;
        oc::BitVector occupy_indicator_flat;

        std::vector<uint32_t> ot_idx;
        std::vector<uint64_t> maskings;
        // Full-layout masks, one per (layer, bin) slot, layer-major:
        // maskings_full[lay * mNumSlots + bin]. Populated only when
        // mFullLayout is set; the equality then runs over every slot.
        std::vector<uint64_t> maskings_full;
        std::vector<seal::Plaintext> ptxts_mask;

        std::vector<std::vector<seal::Plaintext>> ptxts_diags;

        bool mRpmt = false;
        // Full-layout ss-PMT: run the equality over the ENTIRE L x H layout
        // rather than only the n_y occupied slots, and never transmit the
        // occupancy. This closes the layout-disclosure leak (the receiver
        // uses a public band hash, so it can predict its own items' slot
        // columns; occupancy would then be a membership oracle -- see the
        // KKLS follow-up note, Topic B). The sequenced layers are randomly
        // permuted so the greedy front-loading (early layers dense, later
        // layers sparse) does not survive as a positional prior either.
        bool mFullLayout = false;
        uint64_t mOTeBatchSize = 1ull << 19;

    public:
        void sequencing(const std::vector<uint32_t>& start_pos_spacing);

        void init(
            uint32_t n, uint32_t nReceiver,
            sspmtParams ssParams, oc::block seed = oc::OneBlock);

        void rpmt_on() {mRpmt = true;};
        void fullLayoutOn() {mFullLayout = true;};
        auto get_ot_idx() {return ot_idx;};
        u32 getNumLayers() {return mNumLayers;};
        // Number of equality instances actually run: n_y in the compact
        // mode, L * H in the full-layout mode.
        u64 getLayoutSize() {return (u64)mNumLayers * mNumSlots;};

        Proto run(
            const std::vector<oc::block> &Y, 
            Socket &chl);

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

        uint64_t mIndicatorStr;

        std::vector<uint32_t> last_layer_per_bin;
        std::vector<std::vector<uint32_t>> mLayerToBins;

        bool mRpmt = false;
        bool mFullLayout = false;
        uint64_t mOTeBatchSize = 1ull << 19;


    public:
        void init(
            uint32_t n, uint32_t nSender,
            sspmtParams ssParams, oc::block seed = oc::ZeroBlock);

        void rpmt_on() {mRpmt = true;};
        void fullLayoutOn() {mFullLayout = true;};

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

        // Receive the decoded ciphertexts when the layer count is already
        // known (read separately). Split out of recv_decoded_chunks so the
        // full-layout mode can learn the layer count first, then receive the
        // ciphertexts on the base channel while generating GMW triples on a
        // forked channel concurrently.
        Proto recv_decoded_body(
            std::vector<seal::Ciphertext> &decoded_in_he,
            uint32_t numLayers,
            Socket &chl);

        void decrypt(
            const std::vector<seal::Ciphertext> &decoded_in_he, 
            std::vector<uint64_t> &dec_results);        
    };
}
