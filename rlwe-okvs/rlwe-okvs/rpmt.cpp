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
    void RpmtSender::init(uint32_t n, uint32_t nReceiver, uint32_t logp, uint64_t numSlots)
    {
        EncryptionParameters parms(scheme_type::bgv);
        mNumSlots = numSlots;
        parms.set_poly_modulus_degree(mNumSlots);

        mN = n;
        mNreceiver = nReceiver;
        mM = (ceil((1.16 * n) / mNumSlots) + 1) * mNumSlots;
        mW = 134; // 1.16 ; security paraemter = 40
        mNumBatch = mM / mNumSlots;
        
        parms.set_coeff_modulus(CoeffModulus::Create(mNumSlots, { 40, 40, 58, 50}));
        mModulus = PlainModulus::Batching(mNumSlots, logp);
        parms.set_plain_modulus(mModulus);
        
        mContext = make_shared<SEALContext>(parms);
        mBatchEncoder = make_unique<BatchEncoder>(*mContext);
        mEvaluator = make_unique<Evaluator>(*mContext);
    };

    Proto RpmtSender::run(
        const std::vector<oc::block> &Y, 
        std::vector<uint32_t> &ot_idx,
        Socket &chl)
    {
        vector<Plaintext> ptxts_diag;
        vector<Plaintext> ptxts_sdiag;
        vector<uint32_t> bin_sizes;
        preprocess(Y, ptxts_diag, ptxts_sdiag, bin_sizes, ot_idx);

        SEALContext context = *mContext;

        string recvstring;
        co_await chl.recvResize(recvstring);        

        stringstream recvstream;
        recvstream.write(recvstring.data(), recvstring.size());
        
        vector<Ciphertext> encoded_in_he(mNumBatch);
        vector<Ciphertext> encoded_in_he_shift(mNumBatch);
        
        for (size_t i = 0; i < mNumBatch; i++) {
            encoded_in_he[i].unsafe_load(context, recvstream);
        }
        for (size_t i = 0; i < mNumBatch; i++) {
            encoded_in_he_shift[i].unsafe_load(context, recvstream);
        }
        setTimePoint("Sender::Recv ctxts & Serialize");
        
        vector<Ciphertext> decoded_in_he;
        encrypted_decode(encoded_in_he, encoded_in_he_shift,
                         ptxts_diag, ptxts_sdiag, decoded_in_he);

        stringstream sendstream;
        for (size_t i = 0; i < decoded_in_he.size(); i++) {
            auto byte = decoded_in_he[i].save(sendstream);
            if (i == 0) cout << "result ctxt in bytes: " << byte << endl;
        }

        co_await chl.send(decoded_in_he.size());        
        co_await chl.send(move(sendstream.str()));
        co_await chl.send(move(bin_sizes));
        setTimePoint("Sender::Serialize & Send");        
    }

    void RpmtSender::preprocess(
        const std::vector<oc::block> &Y,
        std::vector<seal::Plaintext> &ptxts_diag,
        std::vector<seal::Plaintext> &ptxts_sdiag,
        std::vector<uint32_t> &bin_sizes,
        std::vector<uint32_t> &ot_idx)
    {
        PrimeFieldOkvs okvs;
        // okvs.setTimer(getTimer());
        okvs.init(Y.size(), mM, mW, mModulus);
        vector<uint64_t> bands_flat(mN*mW);
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

        // Optimized Sequencing
        bin_sizes.resize(mNumSlots);
        vector<vector<uint64_t>> T_diag;
        vector<vector<uint64_t>> T_sdiag;
        vector<vector<uint32_t>> item_matrix_mapping;
        vector<uint64_t> zeroRow(mM, 0);
        for(size_t i = 0; i < mN; i++) {
            uint64_t pos = start_pos_spacing[i];
            uint64_t binidx = pos % mNumSlots;
            auto matrix_idx = bin_sizes[binidx];
            bin_sizes[binidx] += 1;
            if (T_diag.size() < matrix_idx + 1) {
                T_diag.push_back(zeroRow);
                T_sdiag.push_back(zeroRow);
                item_matrix_mapping.push_back(vector<uint32_t>(mNumSlots, -1));
            }
            item_matrix_mapping[matrix_idx][binidx] = i;
            bool wrap = false;
            for(size_t w = 0; w < mW; w++){
                if(pos >= mM) {
                    pos = pos - mM;
                    wrap = true;
                }
                if (!wrap) T_diag[matrix_idx][pos] = bands_flat[i*mW+w];
                else T_sdiag[matrix_idx][pos] = bands_flat[i*mW+w];
                pos += mNumSlots;
            }
        }
        auto L = T_diag.size();

        ot_idx.resize(mN);
        auto idx = 0;
        for (size_t i = 0; i < L; i++) {
            for (size_t j = 0; j < mNumSlots; j++) {
                if (item_matrix_mapping[i][j] != -1) {
                    ot_idx[idx++] = item_matrix_mapping[i][j];
                }
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
        for (size_t i = 0; i < decoded_in_he.size(); i++) {
            mEvaluator->mod_switch_to_next_inplace(decoded_in_he[i]);
        }
        setTimePoint("Sender::Encrypted OKVS Decoding");
    }

    void RpmtReceiver::init(uint32_t n, uint32_t nSender, uint32_t logp, uint64_t numSlots)
    {
        EncryptionParameters parms(scheme_type::bgv);
        mNumSlots = numSlots;
        parms.set_poly_modulus_degree(mNumSlots);

        mN = n;
        mNsender = nSender;
        mM = (ceil((1.16 * n) / mNumSlots) + 1) * mNumSlots;
        mW = 134; 
        mNumBatch = mM / mNumSlots;
        mPrng.SetSeed(oc::toBlock(0xDEADBEEF12345678ull, 0xBADC0FFEE0DDF00Dull));
        
        parms.set_coeff_modulus(CoeffModulus::Create(mNumSlots, { 40, 40, 58, 50}));
        mModulus = PlainModulus::Batching(mNumSlots, logp);
        parms.set_plain_modulus(mModulus);
        
        mContext = make_shared<SEALContext>(parms);
        KeyGenerator keygen(*mContext);
        SecretKey secret_key = keygen.secret_key();
        PublicKey public_key;
        keygen.create_public_key(public_key);
        
        mEncryptor = make_unique<Encryptor>(*mContext, public_key);
        mEncryptor->set_secret_key(secret_key);

        mBatchEncoder = make_unique<BatchEncoder>(*mContext);
        mDecryptor = make_unique<Decryptor>(*mContext, secret_key);
    };

    Proto RpmtReceiver::run(
        const std::vector<oc::block> &X, 
        oc::BitVector &results,
        Socket &chl)
    {
        stringstream sendstream;
        encode_and_encrypt(X, sendstream);        
        co_await chl.send(move(sendstream.str()));
        setTimePoint("Receiver::Serialize & Send");
        
        size_t decoded_he_size;
        vector<uint32_t> bin_sizes;
        string recvstring;
        co_await chl.recv(decoded_he_size);
        co_await chl.recvResize(recvstring);
        co_await chl.recvResize(bin_sizes);
        stringstream recvstream;
        recvstream.write(recvstring.data(), recvstring.size());

        setTimePoint("Receiver::Recv back and Serialize");
                
        vector<Ciphertext> decoded_in_he(decoded_he_size);
        SEALContext context = *mContext;
        for (size_t i = 0; i < decoded_he_size; i++) {
            decoded_in_he[i].unsafe_load(context, recvstream);
        }

        decrypt(decoded_in_he, bin_sizes, results);
    }

    void RpmtReceiver::encode_and_encrypt_noserialize(
        const std::vector<oc::block> &X, 
        std::vector<Ciphertext> &ctxts,
        std::vector<Ciphertext> &ctxts_shift
    )
    {
        mIndicatorStr = seal::util::barrett_reduce_64(mPrng.get<uint64_t>(), mModulus);
        vector<uint64_t> val(mN, mIndicatorStr);

        PrimeFieldOkvs okvs;
        // okvs.setTimer(getTimer());
        okvs.init(X.size(), mM, mW, mModulus);

        vector<uint64_t> encoded(mM);

        if (!okvs.encode(X, val, encoded)) throw RTE_LOC;        

        vector<uint64_t> encoded_spacing(mM);
        for (uint32_t i = 0; i < encoded.size(); i++) {
            uint32_t q = i / mNumSlots;
            uint32_t r = i % mNumSlots;
            encoded_spacing[i] = encoded[r * mNumBatch + q];
        }
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
            mEncryptor->encrypt_symmetric(ptxt, ctxts[i]);
        }

        for (uint32_t i = 0; i < mNumBatch - 1; i++) {
            for (uint64_t j = 0 ; j < mNumSlots; j++) {
                plainVec[j] = encoded[i*mNumSlots + j + 1];
            }
            mBatchEncoder->encode(plainVec, ptxt);
            mEncryptor->encrypt_symmetric(ptxt, ctxts_shift[i]);
        }
        for (uint64_t j = 0 ; j < mNumSlots-1; j++) {
            plainVec[j] = encoded[(mNumBatch-1)*mNumSlots + j + 1];
        }
        plainVec[mNumSlots-1] = encoded[0];
        mBatchEncoder->encode(plainVec, ptxt);
        mEncryptor->encrypt_symmetric(ptxt, ctxts_shift[mNumBatch - 1]);
        setTimePoint("Receiver::Encryption");
    }

    void RpmtReceiver::encode_and_encrypt(
            const std::vector<oc::block> &X, 
            stringstream &ctxtstream)
    {
        mIndicatorStr = seal::util::barrett_reduce_64(mPrng.get<uint64_t>(), mModulus);
        vector<uint64_t> val(mN, mIndicatorStr);

        PrimeFieldOkvs okvs;
        // okvs.setTimer(getTimer());
        okvs.init(X.size(), mM, mW, mModulus);

        vector<uint64_t> encoded(mM);

        if (!okvs.encode(X, val, encoded)) throw RTE_LOC;        

        vector<uint64_t> encoded_spacing(mM);
        for (uint32_t i = 0; i < encoded.size(); i++) {
            uint32_t q = i / mNumSlots;
            uint32_t r = i % mNumSlots;
            encoded_spacing[i] = encoded[r * mNumBatch + q];
        }
        encoded = encoded_spacing;

        setTimePoint("Receiver::OKVS Encoding");

        vector<uint64_t> plainVec(mNumSlots);
        Plaintext ptxt;
        
        // TODO: Maybe accelerated with vector iterators
        for (uint32_t i = 0; i < mNumBatch; i++) {
            for (uint64_t j = 0 ; j < mNumSlots; j++) {
                plainVec[j] = encoded[i*mNumSlots + j];
            }
            mBatchEncoder->encode(plainVec, ptxt);
            auto byte = mEncryptor->encrypt_symmetric(ptxt).save(ctxtstream);
            if (i == 0) cout << "encryption in byte: " << byte << endl;
        }

        for (uint32_t i = 0; i < mNumBatch - 1; i++) {
            for (uint64_t j = 0 ; j < mNumSlots; j++) {
                plainVec[j] = encoded[i*mNumSlots + j + 1];
            }
            mBatchEncoder->encode(plainVec, ptxt);
            mEncryptor->encrypt_symmetric(ptxt).save(ctxtstream);
        }
        for (uint64_t j = 0 ; j < mNumSlots-1; j++) {
            plainVec[j] = encoded[(mNumBatch-1)*mNumSlots + j + 1];
        }
        plainVec[mNumSlots-1] = encoded[0];
        mBatchEncoder->encode(plainVec, ptxt);
        mEncryptor->encrypt_symmetric(ptxt).save(ctxtstream);
        setTimePoint("Receiver::Encryption");
    }

    void RpmtReceiver::decrypt(
        const std::vector<seal::Ciphertext> &decoded_in_he, 
        const std::vector<uint32_t> &bin_sizes,
        oc::BitVector &results)
    {
        const size_t L = decoded_in_he.size();
        vector<Plaintext> ptxts(L);
   
        // cout << "Noise Budget: " 
        // << mDecryptor->invariant_noise_budget(decoded_in_he[0]) << endl;

        for(size_t i = 0; i < L; i++){
            mDecryptor->decrypt(decoded_in_he[i], ptxts[i]);
        }

        results.resize(mNsender);
        vector<vector<uint64_t>> decodeVec(L, vector<uint64_t>(mNumSlots, 0));

        auto idx = 0;
        for(size_t i = 0; i < L; i++){
            mBatchEncoder->decode(ptxts[i], decodeVec[i]);
            for(size_t j = 0; j < mNumSlots; j++) {
                if (i < bin_sizes[j]) {
                    results[idx++] = (decodeVec[i][j] == mIndicatorStr) ? 1 : 0;
                }
            }
        }
        setTimePoint("Receiver::Decrypt");
    }
            
}