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
    
    // TODO: Predefine sspmt params in another header,
    // choose among them according to input size only.
    // struct sspmtParams {
    //     u32 heNumSlots = 1 << 13;
    //     std::vector<int> heCoeffModulus = {40, 40, 58, 50};
    //     u32 hePlainModulusBits = 60;
    //     u32 bandWidth = 134;
    //     double bandExpansion = 1.16;
    // };

    struct sspmtParams {
        u32 heNumSlots = 1 << 13;
        std::vector<int> heCoeffModulus;
        u32 hePlainModulusBits;
        u32 bandWidth;
        double bandExpansion;
        void initialize(int n){
            switch(n){
                case (1ull << 16):
                    heCoeffModulus = {40, 34, 55, 50};
                    hePlainModulusBits = 56;
                    bandWidth = 21;
                    bandExpansion = 2.9;
                    break;
                case (1ull << 18):
                    heCoeffModulus = {40, 36, 57, 50};
                    hePlainModulusBits = 58;
                    // bandWidth = 73;
                    // bandExpansion = 1.3;
                    bandWidth = 52;
                    bandExpansion = 1.5;
                    break;
                case (1ull << 20):
                    heCoeffModulus = {40, 38, 59, 50};
                    hePlainModulusBits = 60;
                    // bandWidth = 142;
                    // bandExpansion = 1.15;
                    bandWidth = 109;
                    bandExpansion = 1.2;
                    break;
                case (1ull << 22):
                    heCoeffModulus = {40, 40, 58, 50};
                    hePlainModulusBits = 60; // SEAL의 한계로 60을 써야 함.
                    bandWidth = 116;
                    bandExpansion = 1.2;
                    break;
                default:
                    heCoeffModulus = {40, 40, 58, 60};
                    hePlainModulusBits = 60;
                    bandWidth = 21;
                    bandExpansion = 2.9;
                    break;
            }
        }
    };

    class SspmtSender: public oc::TimerAdapter
    {
        
    public:
        oc::PRNG mPrng;
        Modulus mModulus;
        uint64_t mNumSlots;
        unique_ptr<Evaluator> mEvaluator;
        unique_ptr<BatchEncoder> mBatchEncoder;
        shared_ptr<SEALContext> mContext;
        
        uint32_t mN, mNreceiver, mM, mW, mNumBatch, mWrap;
        uint32_t mW_seq;
        uint32_t mNumLayers;
        std::vector<uint32_t> mItemToLayerIdx;
        std::vector<uint32_t> mItemToBlockIdx;
        std::vector<std::vector<uint32_t>> mLayerBins;
        std::vector<uint32_t> mLayerMinBlock;
        std::vector<uint32_t> mLayerMaxBlock;
        
        std::vector<uint32_t> bin_sizes;
        oc::BitVector occupy_indicator_flat;

        std::vector<uint32_t> ot_idx;
        std::vector<uint64_t> maskings;
        std::vector<seal::Plaintext> ptxts_mask;
        
        std::vector<std::vector<seal::Plaintext>> ptxts_diags;

        bool mSeqOpti = false;
        bool mRpmt = false;
        uint64_t mOTeBatchSize = 1ull << 19; 
        
    // public:
        void sequencing_with_span(
            const std::vector<uint32_t>& start_pos_spacing,
            uint32_t span_blocks          
        );

        void sequencing_naive(
            const std::vector<uint32_t>& start_pos_spacing
        );

        void init(
            uint32_t n, uint32_t nReceiver, 
            sspmtParams ssParams, oc::block seed = oc::OneBlock);

        void rpmt_on() {mRpmt = true;};
        void seqopti_on() {mSeqOpti = true;};
        auto get_ot_idx() {return ot_idx;};

        Proto run(
            const std::vector<oc::block> &Y, 
            Socket &chl);

        Proto run(
            const std::vector<oc::block> &Y, 
            oc::BitVector &results,
            Socket &chl);

        void preprocess(
            const std::vector<oc::block> &Y);

        void encrypted_decode(
            const std::vector<std::vector<seal::Ciphertext>> &encoded_in_he,
            std::vector<seal::Ciphertext> &decoded_in_he);
    };

    class SspmtReceiver: public oc::TimerAdapter
    {

        
    public:
        oc::PRNG mPrng;
        Modulus mModulus;
        uint64_t mNumSlots;
        shared_ptr<SEALContext> mContext;
        unique_ptr<Encryptor> mEncryptor;
        unique_ptr<BatchEncoder> mBatchEncoder;
        unique_ptr<Decryptor> mDecryptor;
        
        uint32_t mN, mNsender, mM, mW, mNumBatch, mWrap;

        uint64_t mIndicatorStr;

        std::vector<uint32_t> bin_sizes;
        std::vector<oc::BitVector> occupy_indicator;

        bool mSeqOpti = false;
        bool mRpmt = false;
        uint64_t mOTeBatchSize = 1ull << 19;
        

    // public:
        void init(
            uint32_t n, uint32_t nSender, 
            sspmtParams ssParams, oc::block seed = oc::ZeroBlock);

        void rpmt_on() {mRpmt = true;};
        void seqopti_on() {mSeqOpti = true;};

        Proto run(
            const std::vector<oc::block> &X, 
            oc::BitVector &results,
            Socket &chl);

        void encode_and_encrypt(
            const std::vector<oc::block> &X, 
            stringstream &ctxtstream);    

        void decrypt(
            const std::vector<seal::Ciphertext> &decoded_in_he, 
            std::vector<uint64_t> &dec_results);        
    };
}