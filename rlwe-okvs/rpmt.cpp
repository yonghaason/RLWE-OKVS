#include "rpmt.h"
#include <memory>
#include "rlwe-okvs/okvs.h"

using namespace std;
using namespace seal;
using namespace oc;

namespace rlweOkvs 
{
    void RpmtSender::init(uint32_t n, uint32_t logp, uint64_t numSlots)
    {
        EncryptionParameters parms(scheme_type::bgv);
        mNumSlots = numSlots;
        parms.set_poly_modulus_degree(mNumSlots);

        mN = n;
        mM = ((ceil(1.16 * n) / mNumSlots) + 1) * mNumSlots;
        mW = 140; // FIXME
        mNumBatch = mM / mNumSlots;
        
        // TODO: Seek the optimal parameter (particularly coeff modulus)
        parms.set_coeff_modulus(CoeffModulus::BFVDefault(mNumSlots));
        mModulus = PlainModulus::Batching(mNumSlots, logp);
        parms.set_plain_modulus(mModulus);
        SEALContext context(parms);

        mBatchEncoder = make_unique<BatchEncoder>(BatchEncoder(context));
        mEvaluator = make_unique<Evaluator>(Evaluator(context));
    };

    Proto RpmtSender::run(const std::vector<oc::block> &Y, Socket &chl)
    {
        vector<Plaintext> ptxts;
        preprocess(Y, ptxts);

        vector<Ciphertext> encoded_in_he(mNumBatch);
        co_await chl.recv(encoded_in_he);

        vector<Ciphertext> decoded_in_he;
        encrypted_decode(encoded_in_he, ptxts, decoded_in_he);
        co_await chl.send(decoded_in_he);
    }

    void RpmtSender::preprocess(
        const std::vector<oc::block> &Y,
        std::vector<seal::Plaintext> &ptxts)
    {
        
    }

    void RpmtSender::encrypted_decode(
        const std::vector<seal::Ciphertext> &encoded_in_he,
        const std::vector<seal::Plaintext> &ptxts,
        std::vector<seal::Ciphertext> &decoded_in_he)
    {

    }

    void RpmtReceiver::init(uint32_t n, uint32_t logp, uint64_t numSlots)
    {
        EncryptionParameters parms(scheme_type::bgv);
        mNumSlots = numSlots;
        parms.set_poly_modulus_degree(mNumSlots);

        mN = n;
        mM = ((ceil(1.16 * n) / mNumSlots) + 1) * mNumSlots;
        mW = 140; // FIXME
        mNumBatch = mM / mNumSlots;
        
        // TODO: Seek the optimal parameter (particularly coeff modulus)
        parms.set_coeff_modulus(CoeffModulus::BFVDefault(mNumSlots));
        mModulus = PlainModulus::Batching(mNumSlots, logp);
        parms.set_plain_modulus(mModulus);
        SEALContext context(parms);

        KeyGenerator keygen(context);
        SecretKey secret_key = keygen.secret_key();

        mEncryptor = make_unique<Encryptor>(Encryptor(context, secret_key));
        mBatchEncoder = make_unique<BatchEncoder>(BatchEncoder(context));
        mDecryptor = make_unique<Decryptor>(Decryptor(context, secret_key));
    };

    Proto RpmtReceiver::run(
        const std::vector<oc::block> &X, 
        oc::BitVector &results,
        Socket &chl)
    {
        vector<Ciphertext> ctxts(mNumBatch);
        encode_and_encrypt(X, ctxts);
        
        // TODO: use seal serialize to compress ctxts (into seed)
        co_await chl.send(std::move(ctxts));

        vector<Ciphertext> decoded_in_he;
        co_await chl.recv(decoded_in_he);

        results.resize(their_size); // initialize 할 때 알려주자
        decrypt(decoded_in_he, results);
    }

    void RpmtReceiver::encode_and_encrypt(
        const std::vector<oc::block> &X, 
        std::vector<seal::Ciphertext> &ctxts
    )
    {
        assert (X.size() == mN);
        mIndicatorStr = mPrng.get<uint64_t>();
        vector<uint64_t> val(mN, mIndicatorStr);
        
        PrimeFieldOkvs okvs;
        okvs.setTimer(getTimer());
        okvs.init(X.size(), mM, mW, mModulus);

        vector<uint64_t> encoded(mM);
        okvs.encode(X, val, encoded);

        vector<uint64_t> plainVec(mNumSlots);
        Plaintext ptxt;
        // TODO: Maybe accelerated with vector iterators
        for (int i = 0; i < mNumBatch; i++) {
            for (int j = 0 ; mNumSlots; j++) {
                plainVec[j] = encoded[i*mNumSlots + j];
            }
            mBatchEncoder->encode(plainVec, ptxt);
            mEncryptor->encrypt(ptxt, ctxts[i]);
        }
    }
}