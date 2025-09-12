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

    class RpmtSender: public oc::TimerAdapter
    {
        oc::PRNG mPrng;
        Modulus mModulus;
        uint64_t mNumSlots;
        unique_ptr<Evaluator> mEvaluator;
        unique_ptr<BatchEncoder> mBatchEncoder;
        unique_ptr<Encryptor> mEncryptor;
        shared_ptr<SEALContext> mContext;
        
        uint32_t mN;
        uint32_t mM;
        uint32_t mW; 
        uint32_t mNumBatch;
        
    public:
        void init(uint32_t n, uint32_t logp, uint64_t numSlots);

        void preprocess(
            const std::vector<oc::block> &Y,
            std::vector<seal::Plaintext> &ptxts_diag,
            std::vector<seal::Plaintext> &ptxts_sdiag);

        void encrypted_decode(
            const std::vector<seal::Ciphertext> &encoded_in_he,
            const std::vector<seal::Ciphertext> &encoded_in_he_shift,
            const std::vector<seal::Plaintext> &ptxts_diag,
            const std::vector<seal::Plaintext> &ptxts_sdiag,
            std::vector<seal::Ciphertext> &decoded_in_he);
        // Proto run(const std::vector<oc::block> &Y, Socket &chl);
    };

    class RpmtReceiver: public oc::TimerAdapter
    {
        oc::PRNG mPrng;
        Modulus mModulus;
        uint64_t mNumSlots;
        shared_ptr<SEALContext> mContext;
        unique_ptr<Encryptor> mEncryptor;
        unique_ptr<BatchEncoder> mBatchEncoder;
        unique_ptr<Decryptor> mDecryptor;
        
        uint32_t mN;
        uint32_t mM;
        uint32_t mW; 
        uint32_t mNumBatch;

        uint64_t mIndicatorStr;
        uint64_t their_size;
        
    public:

        void init(uint32_t n, uint32_t logp, uint64_t numSlots);

        // Proto run(
        //     const std::vector<oc::block> &X, 
        //     oc::BitVector &results,
        //     Socket &chl);

        void encode_and_encrypt(
            const std::vector<oc::block> &X, 
            std::vector<seal::Ciphertext> &ctxts,
            std::vector<seal::Ciphertext> &ctxts_shift
        );

        void decrypt(
            std::vector<seal::Ciphertext> &decoded_in_he, 
            oc::BitVector &results);

        
    };
}