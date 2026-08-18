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

        u32 span_blocks = 0;
        double bandExpansion;

        u32 layerBudget = 0;

        u32 layerBudgetLambda = 40;

        u32 floodBits = 54;

        u32 resolveLayerBudget(u64 n) const;

        u32 resolveSpanBlocks(u64 n) const;

        double spanCostRatio = 1.8e-3;

        void initialize(int n) {
            switch(n){
                case (1ull << 16):
                    floodBits = 50;
                    heCoeffModulus = {50, 58, 60, 50};
                    hePlainModulusBits = 56;
                    bandWidth = 29;
                    bandExpansion = 2.1;
                    break;
                case (1ull << 18):
                    floodBits = 52;
                    heCoeffModulus = {54, 58, 60, 46};
                    hePlainModulusBits = 58;
                    bandWidth = 43;
                    bandExpansion = 1.7;
                    break;
                case (1ull << 20):
                    floodBits = 54;
                    heCoeffModulus = {58, 58, 60, 42};
                    hePlainModulusBits = 60;
                    bandWidth = 44;
                    bandExpansion = 1.7;
                    break;
                case (1ull << 22):
                    floodBits = 56;
                    heCoeffModulus = {60, 60, 60, 38};
                    hePlainModulusBits = 60;
                    bandWidth = 46;
                    bandExpansion = 1.7;
                    break;
                default:
                    heCoeffModulus = {58, 58, 60, 42};
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
        uint32_t mNumLayers;
        uint32_t mNumRealLayers;
        uint32_t mLayerBudget;
        uint32_t mSpanBlocks;
        std::vector<uint32_t> mItemToLayerIdx;
        std::vector<uint32_t> mItemToBlockIdx;
        std::vector<uint32_t> mLayerMinBlock;
        std::vector<uint32_t> mLayerMaxBlock;
        std::vector<uint8_t> mLayerIsPadding;

        std::vector<uint32_t> mSlotToItem;

        std::vector<uint64_t> mMasks;

        std::vector<seal::Ciphertext> mFloodedMasks;

        std::vector<std::vector<seal::Plaintext>> ptxts_diags;

        volePSI::Gmw mGmw;
        bool mSetupDone = false;
        uint64_t mOTeBatchSize = 1ull << 19;

        unique_ptr<Encryptor> mEncryptor;
        seal::parms_id_type mReturnParms;
        u32 mFloodBits;

        void initGmw();

        void addFloodingNoise(seal::Ciphertext& ct);
        void buildFloodedMasks();

    public:

        static uint32_t sequenceLayers(
            const std::vector<uint32_t>& itemBin,
            const std::vector<uint32_t>& itemBlock,
            uint32_t numSlots, uint32_t spanBlocks,
            std::vector<uint32_t>& itemToLayer,
            std::vector<uint32_t>& layerMinBlock,
            std::vector<uint32_t>& layerMaxBlock);

        void sequencing(const std::vector<uint32_t>& start_pos_spacing);

        void init(
            uint32_t n, uint32_t nReceiver,
            sspmtParams ssParams, oc::block seed = oc::OneBlock);

        Proto setup(Socket& chl);

        u32 getNumLayers() {return mNumLayers;};

        u64 getLayoutSize() {return (u64)mNumLayers * mNumSlots;};

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

        PublicKey mPublicKey;

        volePSI::Gmw mGmw;
        bool mSetupDone = false;
        uint64_t mOTeBatchSize = 1ull << 19;

        void initGmw();

    public:
        void init(
            uint32_t n, uint32_t nSender,
            sspmtParams ssParams, oc::block seed = oc::ZeroBlock);

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
