// #pragma once
// #include "cryptoTools/Common/Defines.h"
// #include "cryptoTools/Common/block.h"
// #include "cryptoTools/Common/BitVector.h"
// #include "cryptoTools/Crypto/PRNG.h"
// #include "cryptoTools/Common/Timer.h"
// #include "cryptoTools/Network/Channel.h"
// #include "coproto/coproto.h"

// #include "seal/seal.h"

// using namespace std;
// using namespace seal;
// using namespace oc;

// namespace rlweOkvs
// {
//     using Proto = coproto::task<>;
//     using Socket = coproto::Socket;
    
//     // TODO: Predefine sspmt params in another header,
//     // choose among them according to input size only.
//     struct sspmtParams {
//         u32 heNumSlots = 1 << 13;
//         std::vector<int> heCoeffModulus = {40, 40, 58, 50};
//         u32 hePlainModulusBits = 60;
//         u32 bandWidth = 134;
//         double bandExpansion = 1.16;        
//     };

//     class SspmtSender: public oc::TimerAdapter
//     {
//         oc::PRNG mPrng;
//         Modulus mModulus;
//         uint64_t mNumSlots;
//         unique_ptr<Evaluator> mEvaluator;
//         unique_ptr<BatchEncoder> mBatchEncoder;
//         shared_ptr<SEALContext> mContext;
        
//         uint32_t mN, mNreceiver, mM, mW, mNumBatch, mWrap;

//         std::vector<uint32_t> ot_idx;
//         std::vector<uint64_t> maskings;
//         std::vector<seal::Plaintext> ptxts_mask;
//         std::vector<uint32_t> bin_sizes;

//         bool mRpmt = false;
//         uint64_t mOTeBatchSize = 1ull << 19;
        
//     public:

//         void init(
//             uint32_t n, uint32_t nReceiver, 
//             sspmtParams ssParams, oc::block seed = oc::OneBlock);

//         void rpmt_on() {mRpmt = true;};
//         auto get_ot_idx() {return ot_idx;};

//         Proto run(
//             const std::vector<oc::block> &Y, 
//             Socket &chl);

//         Proto run(
//             const std::vector<oc::block> &Y, 
//             oc::BitVector &results,
//             Socket &chl);

//         void preprocess(
//             const std::vector<oc::block> &Y,
//             std::vector<seal::Plaintext> &ptxts_diag,
//             std::vector<seal::Plaintext> &ptxts_sdiag);

//         void encrypted_decode(
//             const std::vector<seal::Ciphertext> &encoded_in_he,
//             const std::vector<seal::Ciphertext> &encoded_in_he_shift,
//             const std::vector<seal::Plaintext> &ptxts_diag,
//             const std::vector<seal::Plaintext> &ptxts_sdiag,
//             std::vector<seal::Ciphertext> &decoded_in_he);
//     };

//     class SspmtReceiver: public oc::TimerAdapter
//     {
//         oc::PRNG mPrng;
//         Modulus mModulus;
//         uint64_t mNumSlots;
//         shared_ptr<SEALContext> mContext;
//         unique_ptr<Encryptor> mEncryptor;
//         unique_ptr<BatchEncoder> mBatchEncoder;
//         unique_ptr<Decryptor> mDecryptor;
        
//         uint32_t mN, mNsender, mM, mW, mNumBatch, mWrap;

//         uint64_t mIndicatorStr;

//         std::vector<uint32_t> bin_sizes;

//         bool mRpmt = false;
//         uint64_t mOTeBatchSize = 1ull << 19;
        
//     public:

//         void init(
//             uint32_t n, uint32_t nSender, 
//             sspmtParams ssParams, oc::block seed = oc::ZeroBlock);

//         void rpmt_on() {mRpmt = true;};

//         Proto run(
//             const std::vector<oc::block> &X, 
//             oc::BitVector &results,
//             Socket &chl);

//         void encode_and_encrypt(
//             const std::vector<oc::block> &X, 
//             stringstream &ctxtstream);    

//         void decrypt(
//             const std::vector<seal::Ciphertext> &decoded_in_he, 
//             std::vector<uint64_t> &dec_results);        
//     };
// }