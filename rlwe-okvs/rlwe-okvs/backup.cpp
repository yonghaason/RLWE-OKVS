#include "sspmt.h"
#include "okvs.h"
#include "seal/util/uintarithsmallmod.h"
#include "seal/util/numth.h"
#include "../GMW/Gmw.h"
#include <memory>

#include <set>

using namespace std;
using namespace seal;
using namespace oc;
using namespace volePSI;

namespace rlweOkvs 
{
    void SspmtSender::init(
        uint32_t n, uint32_t nReceiver, 
        sspmtParams ssParams, oc::block seed)
    {
        EncryptionParameters parms(scheme_type::bgv);
        mNumSlots = ssParams.heNumSlots;
        parms.set_poly_modulus_degree(mNumSlots);

        mN = n;
        mNreceiver = nReceiver;
        mM = roundUpTo(ssParams.bandExpansion * n, mNumSlots);
        mW = ssParams.bandWidth; // 134 for 1.16
        mNumBatch = mM / mNumSlots;
        mWrap = divCeil(mW * mNumSlots, mM) + 1;
        mPrng.SetSeed(seed);        
        
        parms.set_coeff_modulus(CoeffModulus::Create(mNumSlots, ssParams.heCoeffModulus));
        mModulus = PlainModulus::Batching(mNumSlots, ssParams.hePlainModulusBits);
        parms.set_plain_modulus(mModulus);
        
        mContext = make_shared<SEALContext>(parms);
        mBatchEncoder = make_unique<BatchEncoder>(*mContext);
        mEvaluator = make_unique<Evaluator>(*mContext);
    };

    Proto SspmtSender::run(
        const std::vector<oc::block> &Y, 
        Socket &chl)
    {
        vector<Plaintext> ptxts_diag;
        vector<Plaintext> ptxts_sdiag;
        preprocess(Y, ptxts_diag, ptxts_sdiag);

        SEALContext context = *mContext;

        string recvstring;
        co_await chl.recvResize(recvstring);
        cout << "Sender receives " << mNumBatch << " X " << mWrap 
            << " Ctxts of OKVS Encoding, "
            << recvstring.size() << " Bytes" << endl;

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
            // if (i == 0) cout << "result ctxt in bytes: " << byte << endl;
        }

        cout << "Sender sends " << decoded_in_he.size() << " decoded ctxts, "
             << sendstream.str().size() << " Bytes" << endl;

        co_await chl.send(decoded_in_he.size());        
        co_await chl.send(move(sendstream.str()));
        co_await chl.send(move(bin_sizes));
        setTimePoint("Sender::Serialize & Send");        
    }

    Proto SspmtSender::run(
        const std::vector<oc::block> &Y, 
        oc::BitVector &results,
        Socket &chl)
    {
        assert(!mRpmt && "Sender obtains result only when ssPMT.");
        co_await run(Y, chl);

        // GMW with maskings
        u64 keyBitLength = 40 + oc::log2ceil(Y.size());
        u64 keyByteLength = oc::divCeil(keyBitLength, 8);

        oc::Matrix<u8> gmwin;
        gmwin.resize(Y.size(), keyByteLength, oc::AllocType::Uninitialized);
        for (size_t i = 0; i < Y.size(); i++) {
            memcpy(&gmwin(i, 0), &maskings[i], keyByteLength);
        }
        
        Gmw gmw;
        gmw.setTimer(getTimer());
        auto cir = isZeroCircuit(keyBitLength);
        gmw.init(Y.size(), cir, 1, mOTeBatchSize, 0, mPrng.get());
        gmw.setInput(0, gmwin);

        co_await gmw.run(chl);
        auto rr = gmw.getOutputView(0);
        results.resize(Y.size());
        std::copy(rr.begin(), rr.end(), results.data());
        setTimePoint("Sender::Online GMW");
    }

    void SspmtSender::preprocess(
        const std::vector<oc::block> &Y,
        std::vector<seal::Plaintext> &ptxts_diag,
        std::vector<seal::Plaintext> &ptxts_sdiag)
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
        std::vector<uint64_t> T_diag_flat;    // L * mM
        std::vector<uint64_t> T_sdiag_flat;   // L * mM
        
        std::vector<std::vector<uint32_t>> bins_items(mNumSlots);
        std::vector<uint32_t> item_binidx(mN);
        std::vector<uint32_t> item_matrix_idx(mN);

        for (uint32_t i = 0; i < mN; ++i) {
            uint32_t pos = start_pos_spacing[i];
            uint32_t binidx = pos % mNumSlots;
            item_binidx[i] = binidx;
            item_matrix_idx[i] = bins_items[binidx].size();
            bins_items[binidx].push_back(i);
        }

        bin_sizes.resize(mNumSlots);
        for (size_t i = 0; i < mNumSlots; i++) {
            bin_sizes[i] = bins_items[i].size();
        }

        size_t L = 0;
        for (const auto& v : bins_items) L = std::max(L, v.size());

        T_diag_flat.assign(L * mM, 0);
        T_sdiag_flat.assign(L * mM, 0);

        auto* diag = T_diag_flat.data();
        auto* sdiag = T_sdiag_flat.data();
        const uint64_t* bands = bands_flat.data();

        for (uint32_t i = 0; i < mN; ++i) {
            const uint32_t pos = start_pos_spacing[i];
            const uint32_t row = item_matrix_idx[i];
            const uint64_t* src = bands + i * mW;

            uint32_t room = (mM - pos + mNumSlots - 1) / mNumSlots;
            size_t take0 = std::min<uint64_t>(mW, room);

            uint32_t p = pos;
            for (size_t w = 0; w < take0; ++w, p += mNumSlots) {
                diag[row * mM + p] = src[w];
            }

            if (take0 < mW) {
                uint32_t p2 = pos + take0 * mNumSlots - mM;
                for (size_t w = take0; w < mW; ++w, p2 += mNumSlots) {
                    sdiag[row * mM + p2] = src[w];
                }
            }
        }
        
        ot_idx.resize(mN);
        size_t idx = 0;

        if (mRpmt) {
            for (size_t mat = 0; mat < L; ++mat) {
                for (uint32_t bin = 0; bin < mNumSlots; ++bin) {
                    const auto& vec = bins_items[bin];
                    if (mat < vec.size())
                        ot_idx[idx++] = vec[mat];
                }
            }
        }
        else {
            maskings.resize(mN);
            ptxts_mask.resize(L);
            vector<vector<uint64_t>> raw_masks(L);            
            for (size_t mat = 0; mat < L; mat++) {
                raw_masks[mat].resize(mNumSlots);
                mPrng.get<uint64_t>(raw_masks[mat]);
                for (uint32_t bin = 0; bin < mNumSlots; bin++) {
                    raw_masks[mat][bin] = seal::util::barrett_reduce_64(raw_masks[mat][bin], mModulus);
                    const auto& vec = bins_items[bin];
                    if (mat < vec.size()) {
                        ot_idx[idx] = vec[mat];
                        maskings[idx++] = raw_masks[mat][bin];
                    }
                }
                mBatchEncoder->encode(raw_masks[mat], ptxts_mask[mat]);
            }
        }

        //Batch  
        vector<uint64_t> plainVec(mNumSlots);
        ptxts_diag.resize(mNumBatch*L);
        ptxts_sdiag.resize(mNumBatch*L);
    
        size_t outIdx = 0;
        for (size_t i = 0; i < L; ++i) {
            const size_t row_base = i * mM;
            for (uint32_t j = 0; j < mNumBatch; ++j) {
                const size_t off = row_base + j * mNumSlots;
                std::copy_n(&T_diag_flat[off], mNumSlots, plainVec.data());
                mBatchEncoder->encode(plainVec, ptxts_diag[outIdx]);
                std::copy_n(&T_sdiag_flat[off], mNumSlots, plainVec.data());
                mBatchEncoder->encode(plainVec, ptxts_sdiag[outIdx++]);
            }
        }
        setTimePoint("Sender::Sequencing & HE Encode");
    }

    void SspmtSender::encrypted_decode(
        const std::vector<seal::Ciphertext> &encoded_in_he,
        const std::vector<seal::Ciphertext> &encoded_in_he_shift,
        const std::vector<seal::Plaintext> &ptxts_diag,
        const std::vector<seal::Plaintext> &ptxts_sdiag,
        std::vector<seal::Ciphertext> &decoded_in_he)
    {
        u32 L = 0;
        for (auto& s: bin_sizes) L = max(L, s);

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

        if (!mRpmt) {
            for (size_t i = 0; i < L; ++i) {
                mEvaluator->add_plain_inplace(decoded_in_he[i], ptxts_mask[i]);
            }
        }

        for (size_t i = 0; i < decoded_in_he.size(); i++) {
            mEvaluator->mod_switch_to_next_inplace(decoded_in_he[i]);
        }
        setTimePoint("Sender::Encrypted OKVS Decoding");
    }

    void SspmtReceiver::init(
        uint32_t n, uint32_t nSender, 
        sspmtParams ssParams, oc::block seed)
    {
        EncryptionParameters parms(scheme_type::bgv);
        mNumSlots = ssParams.heNumSlots;
        parms.set_poly_modulus_degree(mNumSlots);

        mN = n;
        mNsender = nSender;
        mM = roundUpTo(ssParams.bandExpansion * n, mNumSlots);
        mW = ssParams.bandWidth; // 134 for 1.16
        mNumBatch = mM / mNumSlots;
        mWrap = divCeil(mW * mNumSlots, mM) + 1;
        mPrng.SetSeed(seed);

        parms.set_coeff_modulus(CoeffModulus::Create(mNumSlots, ssParams.heCoeffModulus));
        mModulus = PlainModulus::Batching(mNumSlots, ssParams.hePlainModulusBits);
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

    Proto SspmtReceiver::run(
        const std::vector<oc::block> &X, 
        oc::BitVector &results,
        Socket &chl)
    {
        stringstream sendstream;
        encode_and_encrypt(X, sendstream);        
        co_await chl.send(move(sendstream.str()));
        setTimePoint("Receiver::Serialize & Send");
        
        size_t decoded_he_size;
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

        vector<uint64_t> dec_results;
        decrypt(decoded_in_he, dec_results);

        if (mRpmt) {
            results.resize(mNsender);
            for(size_t i = 0; i < mNsender; i++){
                results[i] = (dec_results[i] == mIndicatorStr) ? 1 : 0;
            }
        }
        else {
            // GMW with dec_results
            u64 keyBitLength = 40 + oc::log2ceil(mNsender);
            u64 keyByteLength = oc::divCeil(keyBitLength, 8);

            oc::Matrix<u8> gmwin;
            
            gmwin.resize(mNsender, keyByteLength, oc::AllocType::Uninitialized);
            for (size_t i = 0; i < mNsender; i++) {
                auto tmp = util::sub_uint_mod(dec_results[i], mIndicatorStr, mModulus);
                memcpy(&gmwin(i, 0), &tmp, keyByteLength);
            }
            
            Gmw gmw;
            gmw.setTimer(getTimer());
            auto cir = isZeroCircuit(keyBitLength);
            gmw.init(mNsender, cir, 1, mOTeBatchSize, 1, mPrng.get());
            gmw.setInput(0, gmwin);

            co_await gmw.run(chl);
            auto rr = gmw.getOutputView(0);
            results.resize(mNsender);
            std::copy(rr.begin(), rr.end(), results.data());
            setTimePoint("Receiver::Online GMW");
        }
    }

    void SspmtReceiver::encode_and_encrypt(
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
            // if (i == 0) cout << "encryption in byte: " << byte << endl;
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

    void SspmtReceiver::decrypt(
        const std::vector<seal::Ciphertext> &decoded_in_he,
        std::vector<uint64_t> &dec_results)
    {
        const size_t L = decoded_in_he.size();
        vector<Plaintext> ptxts(L);
   
        // cout << "Noise Budget: " 
        // << mDecryptor->invariant_noise_budget(decoded_in_he[0]) << endl;

        for(size_t i = 0; i < L; i++){
            mDecryptor->decrypt(decoded_in_he[i], ptxts[i]);
        }

        dec_results.resize(mNsender);
        vector<uint64_t> decodeVec(mNumSlots);

        auto idx = 0;
        for(size_t i = 0; i < L; i++){
            mBatchEncoder->decode(ptxts[i], decodeVec);
            for(size_t j = 0; j < mNumSlots; j++) {
                if (i < bin_sizes[j]) {
                    dec_results[idx++] = decodeVec[j];
                }
            }
        }
        setTimePoint("Receiver::Decrypt");
    }            
}