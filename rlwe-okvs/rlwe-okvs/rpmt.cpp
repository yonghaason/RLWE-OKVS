#include "rpmt.h"
#include "okvs.h"
#include "seal/util/uintarithsmallmod.h"
#include "seal/util/numth.h"
#include <memory>

#include <set>

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
        mM = (ceil((1.16 * n) / mNumSlots) + 1) * mNumSlots;
        mW = 134; // 1.16 ; security paraemter = 40
        mNumBatch = mM / mNumSlots;
        
        // TODO: Seek the optimal parameter (particularly coeff modulus)
        // parms.set_coeff_modulus(CoeffModulus::BFVDefault(mNumSlots));
        parms.set_coeff_modulus(CoeffModulus::Create(mNumSlots, { 48, 48, 48, 50}));
        mModulus = PlainModulus::Batching(mNumSlots, logp);
        parms.set_plain_modulus(mModulus);
        
        mContext = make_shared<SEALContext>(parms);
        mBatchEncoder = make_unique<BatchEncoder>(*mContext);
        mEvaluator = make_unique<Evaluator>(*mContext);
    };

    // Proto RpmtSender::run(const std::vector<oc::block> &Y, Socket &chl)
    // {
    //     vector<Plaintext> ptxts;
    //     preprocess(Y, ptxts);

    //     vector<Ciphertext> encoded_in_he(mNumBatch);
    //     co_await chl.recv(encoded_in_he);

    //     vector<Ciphertext> decoded_in_he;
    //     encrypted_decode(encoded_in_he, ptxts, decoded_in_he);
    //     co_await chl.send(decoded_in_he);
    // }

    void RpmtSender::preprocess(
        const std::vector<oc::block> &Y,
        std::vector<seal::Plaintext> &ptxts_diag,
        std::vector<seal::Plaintext> &ptxts_sdiag)
    {
        //Make matrix Y
        assert (Y.size() == mN);

        PrimeFieldOkvs okvs;
        // okvs.setTimer(getTimer());
        okvs.init(Y.size(), mM, mW, mModulus);
        vector<uint64_t> bands_flat(mN*mW); // this is matrix Y
        vector<uint32_t> start_pos(mN);
        okvs.generate_band(Y, bands_flat, start_pos, oc::ZeroBlock);

        setTimePoint("Sender::Generate Bands");

        vector<uint32_t> start_pos_spacing(mN);
        for (uint32_t i = 0; i < mN; i++) {
            auto position = start_pos[i];
            uint32_t q = position / mNumBatch;
            uint32_t r = position % mNumBatch;
            start_pos_spacing[i] = r * mNumSlots + q;
        }

        //get number of t_s set
        vector<uint32_t> start_freq(mNumSlots, 0);
        uint64_t L = 0; // number of t_s set

        uint32_t* f = start_freq.data();

        for(size_t i = 0; i < mN; ++i){
            uint64_t r = start_pos_spacing[i] % mNumSlots; 
            uint32_t c = ++f[static_cast<size_t>(r)];
            if(c > L) L = c;
        }

        //Sequencing
        vector<vector<uint64_t>> T_diag(L, vector<uint64_t>(mM, 0));
        vector<vector<uint64_t>> T_sdiag(L, vector<uint64_t>(mM, 0));
        vector<oc::BitVector> filled(L);
        for (size_t i = 0; i < L; i++) {
            filled[i].resize(mNumSlots);
        }
        
        // TODO: Fast Sequencing? Filled 쓰지말고 current ball num 기억하면 될듯
        for(size_t i = 0; i < mN; i++){
            uint64_t binidx = start_pos_spacing[i] % mNumSlots;
            uint64_t ballnum = 0;
            while(true){
                if(!filled[ballnum][binidx]){
                    uint32_t pos = start_pos_spacing[i];
                    bool wrap = false;
                    // TODO: If 없이 optimize 가능 (pos -> mM 얼마나 남은지 계산해서)
                    for(size_t w = 0; w < mW; w++){
                        if(pos >= mM) {
                            pos = pos - mM;
                            wrap = true;
                        }
                        if (wrap) T_sdiag[ballnum][pos] = bands_flat[i*mW+w];
                        else T_diag[ballnum][pos] = bands_flat[i*mW+w];
                        pos += mNumSlots;
                    }
                    filled[ballnum][binidx] = 1;
                    break;
                }
                else ballnum++;
            }
        }
        setTimePoint("Sender::Sequencing");

        //Batch  
        vector<uint64_t> plainVec(mNumSlots);
        ptxts_diag.resize(mNumBatch*L);
        ptxts_sdiag.resize(mNumBatch*L);
        
        size_t outIdx = 0;

        for(uint64_t i = 0; i < L; i++){
            for(uint32_t j = 0 ; j < mNumBatch; j++){
                const size_t offset = j * mNumSlots;
                for(uint64_t k = 0; k < mNumSlots; k++){
                    plainVec[k] = T_diag[i][offset+k];
                }
                mBatchEncoder->encode(plainVec, ptxts_diag[outIdx]);
                for(uint64_t k = 0; k < mNumSlots; k++){
                    plainVec[k] = T_sdiag[i][offset+k];
                }
                mBatchEncoder->encode(plainVec, ptxts_sdiag[outIdx++]);
            }
        }
        setTimePoint("Sender::BatchEncode");
    }

    void RpmtSender::encrypted_decode(
        const std::vector<seal::Ciphertext> &encoded_in_he,
        const std::vector<seal::Ciphertext> &encoded_in_he_shift,
        const std::vector<seal::Plaintext> &ptxts_diag,
        const std::vector<seal::Plaintext> &ptxts_sdiag,
        std::vector<seal::Ciphertext> &decoded_in_he)
    {
        //PlainMult, CtxtAdd
        const size_t L = ptxts_diag.size() / mNumBatch;

        decoded_in_he.resize(L);

        for (size_t i = 0; i < L; ++i) {
            bool initialized = false;
            Ciphertext tmp;

            for (size_t j = 0; j < mNumBatch; ++j) {
                const size_t idx = i * mNumBatch + j;
                const Plaintext &pt_diag = ptxts_diag[idx];
                const Plaintext &pt_sdiag = ptxts_sdiag[idx];

                if (!pt_diag.is_zero()) {
                    if (!initialized) {
                        mEvaluator->multiply_plain(encoded_in_he[j], pt_diag, decoded_in_he[i]);
                        initialized = true;
                    }
                    else {    
                        mEvaluator->multiply_plain(encoded_in_he[j], pt_diag, tmp);
                        mEvaluator->add_inplace(decoded_in_he[i], tmp);
                    }
                }

                if (!pt_sdiag.is_zero()) {
                    if (!initialized) {
                        mEvaluator->multiply_plain(encoded_in_he_shift[j], pt_sdiag, decoded_in_he[i]);
                        initialized = true;
                    }
                    else {    
                        mEvaluator->multiply_plain(encoded_in_he_shift[j], pt_sdiag, tmp);
                        mEvaluator->add_inplace(decoded_in_he[i], tmp);
                    }
                }
            }
        }
        
        setTimePoint("Sender::Encrypted OKVS Decoding");
    }

    void RpmtReceiver::init(uint32_t n, uint32_t logp, uint64_t numSlots)
    {
        EncryptionParameters parms(scheme_type::bgv);
        mNumSlots = numSlots;
        parms.set_poly_modulus_degree(mNumSlots);

        mN = n;
        mM = (ceil((1.16 * n) / mNumSlots) + 1) * mNumSlots;
        mW = 134; // FIXME
        mNumBatch = mM / mNumSlots;
        mPrng.SetSeed(oc::toBlock(0xDEADBEEF12345678ull, 0xBADC0FFEE0DDF00Dull));
        
        // TODO: Seek the optimal parameter (particularly coeff modulus)
        // parms.set_coeff_modulus(CoeffModulus::BFVDefault(mNumSlots));
        parms.set_coeff_modulus(CoeffModulus::Create(mNumSlots, { 48, 48, 48, 50}));
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
        std::vector<seal::Ciphertext> &ctxts,
        std::vector<seal::Ciphertext> &ctxts_shift
    )
    {
        assert (X.size() == mN);
        mIndicatorStr = seal::util::barrett_reduce_64(mPrng.get<uint64_t>(), mModulus);
        vector<uint64_t> val(mN, mIndicatorStr);

        PrimeFieldOkvs okvs;
        // okvs.setTimer(getTimer());
        okvs.init(X.size(), mM, mW, mModulus);

        vector<uint64_t> encoded(mM);

        if (!okvs.encode(X, val, encoded)) throw RTE_LOC;
        
        // vector<uint64_t> decoded(mN);
        // okvs.decode(X, encoded, decoded);

        // for (size_t i = 0; i < val.size(); i++) {
        //     if (val[i] != decoded[i]) {
        //         throw RTE_LOC;
        //     }
        // }
        // std::cout << "Encode is correct" << std::endl;

        vector<uint64_t> encoded_spacing(mM);
        for (uint32_t i = 0; i < encoded.size(); i++) {
            uint32_t q = i / mNumSlots;
            uint32_t r = i % mNumSlots;
            encoded_spacing[i] = encoded[r * mNumBatch + q];
        }

        // vector<uint64_t> bands_flat(mN*mW);
        // vector<uint32_t> start_pos(mN);
        // okvs.generate_band(X, bands_flat, start_pos, oc::ZeroBlock);
        // vector<uint32_t> start_pos_spacing(mN);
        // for (uint32_t i = 0; i< start_pos.size(); i++) {
        //     auto position = start_pos[i];
        //     uint32_t q = position / mNumBatch;
        //     uint32_t r = position % mNumBatch;
        //     start_pos_spacing[i] = r * mNumSlots + q;
        // }
        // mat_vec_mult_spaced_band_flat(bands_flat, start_pos_spacing, mNumSlots, mModulus, encoded_spacing, decoded);

        // auto cnt = 0;
        // for (size_t i = 0; i < val.size(); i++) {
        //     if (val[i] != decoded[i]) {
        //         cout << i << ": " << val[i] << ", " << decoded[i] << endl;
        //         const uint64_t *__restrict row = bands_flat.data() + static_cast<size_t>(i) * mW;

        //         auto position = start_pos_spacing[i];

        //         for (uint32_t j = 0; j < mW; ++j){   
        //             if (position >= mM) position = position - mM + 1;
        //             if (encoded[start_pos[i] + j] != encoded_spacing[position]) {
        //                 std::cout << encoded[start_pos[i] + j] << "(before spacing) v.s. " << encoded_spacing[position] << "(after spacing)" << std::endl;
        //             }
        //             position += mNumSlots;
        //         }
        //         cnt++;
        //     }
        // }
        // if (cnt) throw RTE_LOC;
        // std::cout << "Spacing is correct" << std::endl;

        encoded = encoded_spacing;

        setTimePoint("Receiver::OKVS Encoding");

        ctxts.resize(static_cast<size_t>(mNumBatch));
        ctxts_shift.resize(static_cast<size_t>(mNumBatch));
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

        // TODO: Fast Shift?
        for (uint32_t i = 0; i < mNumBatch - 1; i++) {
            for (uint64_t j = 0 ; j < mNumSlots; j++) {
                plainVec[j] = encoded[i*mNumSlots + j + 1];
            }
            mBatchEncoder->encode(plainVec, ptxt);
            mEncryptor->encrypt(ptxt, ctxts_shift[i]);
        }
        for (uint64_t j = 0 ; j < mNumSlots-1; j++) {
            plainVec[j] = encoded[(mNumBatch-1)*mNumSlots + j + 1];
        }
        plainVec[mNumSlots-1] = encoded[0];
        mBatchEncoder->encode(plainVec, ptxt);
        mEncryptor->encrypt(ptxt, ctxts_shift[mNumBatch - 1]);
        setTimePoint("Receiver::Encryption");
    }

    void RpmtReceiver::decrypt(
        std::vector<seal::Ciphertext> &decoded_in_he, 
        oc::BitVector &results)
    {
        //Decryption
        const size_t L = decoded_in_he.size();
        vector<Plaintext> ptxts(L);

        
        cout << "Noise Budget: "
        << mDecryptor->invariant_noise_budget(decoded_in_he[0]) << endl;

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
        setTimePoint("Receiver::Decrypt");
    }
            
}