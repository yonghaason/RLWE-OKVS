#include "rpmt.h"
#include "okvs.h"
#include <memory>


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
        mW = 134; // 1.16 ; security paraemter = 40
        mNumBatch = mM / mNumSlots;
        
        // TODO: Seek the optimal parameter (particularly coeff modulus)
        parms.set_coeff_modulus(CoeffModulus::BFVDefault(mNumSlots));
        mModulus = PlainModulus::Batching(mNumSlots, logp);
        parms.set_plain_modulus(mModulus);
        
        mContext = make_shared<SEALContext>(parms);
        mBatchEncoder = make_unique<BatchEncoder>(*mContext);
        mEvaluator = make_unique<Evaluator>(*mContext);
    };

    // Proto RpmtSender::run(const std::vector<oc::block> &Y, Socket &chl)
    // {
    //     vector<Plaintext> ptxts;=
    //     preprocess(Y, ptxts);

    //     vector<Ciphertext> encoded_in_he(mNumBatch);
    //     co_await chl.recv(encoded_in_he);

    //     vector<Ciphertext> decoded_in_he;
    //     encrypted_decode(encoded_in_he, ptxts, decoded_in_he);
    //     co_await chl.send(decoded_in_he);
    // }

    /* 구현해야 할 부분 (8.21) */
    void RpmtSender::preprocess(
        const std::vector<oc::block> &Y,
        std::vector<seal::Plaintext> &ptxts)
    {
        //Make matrix Y
        assert (Y.size() == mN);

        PrimeFieldOkvs okvs;
        okvs.setTimer(getTimer());
        okvs.init(Y.size(), mM, mW, mModulus);
        vector<uint64_t> bands_flat(mN*mW); // this is matrix Y
        vector<uint32_t> start_pos(mN);
        okvs.generate_band(Y, bands_flat, start_pos, oc::ZeroBlock);

        //get number of t_s set
        vector<uint32_t> start_freq(mNumSlots, 0);
        uint64_t L = 0; // number of t_s set

        uint32_t* f = start_freq.data();

        for(size_t i = 0; i < mN; ++i){
            uint64_t r = start_pos[i] % mNumSlots; 
            //start_pos[i] = r;
            uint32_t c = ++f[static_cast<size_t>(r)];
            if(c > L) L = c;
        }

        //Sequencing
        setTimePoint("Start sequencing");
        vector<vector<uint64_t>> T_s(L, vector<uint64_t>(mM, 0));
        oc::BitVector filled(L*mNumSlots);

        for(size_t i = 0; i < mN; i++){
            uint64_t r = start_pos[i] % mNumSlots;
            uint64_t j = 0;
            while(true){
                if(!filled[j*mNumSlots+r]){
                    for(size_t w = 0; w < mW; w++){
                        uint32_t p = start_pos[i] + w*mNumSlots;
                        if(p >= mM) p -= mM;
                        T_s[j][p] = bands_flat[i*mW+w];
                    }
                    filled[j*mNumSlots+r] = 1;
                    break;
                }
                else j++;
            }
        }
        setTimePoint("End sequencing");
        
        //Batch  
        vector<uint64_t> plainVec(mNumSlots);
        ptxts.resize(mNumBatch*L);
        size_t outIdx = 0;

        for(uint64_t i = 0; i < L; i++){
            for(uint32_t j = 0 ; j < mNumBatch; j++){
                const size_t offset = j * mNumSlots;
                for(uint64_t k = 0; k < mNumSlots; k++){
                    plainVec[k] = T_s[i][offset+k];
                }
                mBatchEncoder->encode(plainVec, ptxts[outIdx++]);
            }
        }
    }

    /* 구현해야 할 부분 (8.21) */
    void RpmtSender::encrypted_decode(
        const std::vector<seal::Ciphertext> &encoded_in_he,
        const std::vector<seal::Plaintext> &ptxts,
        std::vector<seal::Ciphertext> &decoded_in_he)
    {
        //PlainMult, CtxtAdd
        const size_t L = ptxts.size() / mNumBatch;
        decoded_in_he.resize(L);

        setTimePoint("start decoding(PlainMult&CtxtAdd)");
        for (size_t i = 0; i < L; ++i) {
            Ciphertext acc, tmp;
            bool any_term = false;

            // 모든 블록에 대해: 0-평문은 곱하지 않고 건너뜀
            for (size_t b = 0; b < mNumBatch; ++b) {
                const size_t idx = i * mNumBatch + b;
                const Plaintext &pt = ptxts[idx];

                if (pt.is_zero()) continue;   // ⬅ 0-평문 skip (중요)

                mEvaluator->multiply_plain(encoded_in_he[b], pt, tmp);

                if (!any_term) {
                    acc = std::move(tmp);
                    any_term = true;
                } else {
                    mEvaluator->add_inplace(acc, tmp);
                }
            }

            if (!any_term) {
                // 이 행의 모든 항이 0 → 결과는 0이어야 함.
                // 하지만 투명 0을 만들면 SEAL이 막으므로, 공개키로 '암호화된 0'을 생성해 반환.
                if (!mEncryptor)
                    throw std::runtime_error("RpmtSender: public key not set; cannot emit encrypted zero");
                mEncryptor->encrypt_zero(acc);                 // ⬅ non-transparent zero
            }

            decoded_in_he[i] = std::move(acc);
        }
        setTimePoint("end decoding(PlainMult&CtxtAdd)");
    }

    void RpmtReceiver::init(uint32_t n, uint32_t logp, uint64_t numSlots)
    {
        EncryptionParameters parms(scheme_type::bgv);
        mNumSlots = numSlots;
        parms.set_poly_modulus_degree(mNumSlots);

        mN = n;
        mM = ((ceil(1.16 * n) / mNumSlots) + 1) * mNumSlots;
        mW = 134; // FIXME
        mNumBatch = mM / mNumSlots;
        mPrng.SetSeed(oc::toBlock(0x1234567890ABCDEFull, 0x0FEDCBA098765432ull));
        
        // TODO: Seek the optimal parameter (particularly coeff modulus)
        parms.set_coeff_modulus(CoeffModulus::BFVDefault(mNumSlots));
        mModulus = PlainModulus::Batching(mNumSlots, logp);
        parms.set_plain_modulus(mModulus);
        SEALContext context(parms);
        KeyGenerator keygen(context);
        SecretKey secret_key = keygen.secret_key();
        PublicKey public_key;
        keygen.create_public_key(public_key);
        
        mEncryptor = make_unique<Encryptor>(context, public_key);
        mBatchEncoder = make_unique<BatchEncoder>(context);
        mDecryptor = make_unique<Decryptor>(context, secret_key);
    };

    // Proto RpmtReceiver::run(
    //     const std::vector<oc::block> &X, 
    //     oc::BitVector &results,
    //     Socket &chl)
    // {
    //     vector<Ciphertext> ctxts(mNumBatch);
    //     encode_and_encrypt(X, ctxts);
        
    //     // TODO: use seal serialize to compress ctxts (into seed)
    //     co_await chl.send(std::move(ctxts));

    //     vector<Ciphertext> decoded_in_he;
    //     co_await chl.recv(decoded_in_he);

    //     results.resize(their_size); // initialize 할 때 알려주자
    //     decrypt(decoded_in_he, results);
    // }

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

        ctxts.resize(static_cast<size_t>(mNumBatch));
        vector<uint64_t> plainVec(mNumSlots);
        Plaintext ptxt;
        
        // TODO: Maybe accelerated with vector iterators
        for (uint32_t i = 0; i < mNumBatch; i++) {
            for (uint64_t j = 0 ; j < mNumSlots; j++) {
                plainVec[j] = encoded[i*mNumSlots + j];
            }
            mBatchEncoder->encode(plainVec, ptxt);
            mEncryptor->encrypt(ptxt, ctxts[i]);
        }
    }

    void RpmtReceiver::decrypt(
        std::vector<seal::Ciphertext> &decoded_in_he, 
        oc::BitVector &results)
    {
        //Decryption
        const size_t L = decoded_in_he.size();
        vector<Plaintext> ptxts(L);

        for(size_t i = 0; i < L; i++){
            mDecryptor->decrypt(decoded_in_he[i], ptxts[i]);
        }

        //Decode(Unbatch) & get results for OT
        results.resize(L*mNumSlots);
        vector<vector<uint64_t>> decodeVec(L, vector<uint64_t>(mNumSlots, 0));

        for(size_t i = 0; i < L; i++){
            mBatchEncoder->decode(ptxts[i], decodeVec[i]);
            for(size_t j = 0; j < mNumSlots; j++){
                results[i*mNumSlots+j] = (decodeVec[i][j] == mIndicatorStr) ? 1 : 0;
            }
        }
    }
            
}