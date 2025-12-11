#include "sspmt.h"

#include <memory>
#include <set>

#include "../GMW/Gmw.h"
#include "okvs.h"
#include "seal/util/numth.h"
#include "seal/util/uintarithsmallmod.h"

using namespace std;
using namespace seal;
using namespace oc;
using namespace volePSI;

namespace rlweOkvs {
void SspmtSender::sequencing(const std::vector<uint32_t> &start_pos_spacing) {
  std::vector<uint32_t> item_binidx(mN);
  uint32_t max_block = 0;

  for (uint32_t i = 0; i < mN; ++i) {
    uint32_t pos = start_pos_spacing[i];
    uint32_t bin = pos % mNumSlots;
    uint32_t blk = pos / mNumSlots;

    item_binidx[i] = bin;
    mItemToBlockIdx[i] = blk;
    if (blk > max_block) max_block = blk;
  }

  std::vector<std::vector<uint32_t>> block_items(max_block + 1);
  for (uint32_t i = 0; i < mN; ++i) {
    block_items[mItemToBlockIdx[i]].push_back(i);
  }

  struct Layer {
    uint32_t min_block;
    uint32_t max_block;
    std::vector<uint8_t> used_bins;
  };
  std::vector<Layer> layers;

  for (uint32_t blk = 0; blk <= max_block; ++blk) {
    auto &bucket = block_items[blk];
    for (uint32_t idx : bucket) {
      uint32_t r = item_binidx[idx];      // bin (residue)
      uint32_t j = mItemToBlockIdx[idx];  // block
      bool placed = false;
      for (uint32_t li = 0; li < layers.size(); ++li) {
        Layer &L = layers[li];

        if (L.used_bins[r]) {
          continue;
        }

        uint32_t new_min = std::min(L.min_block, j);
        uint32_t new_max = std::max(L.max_block, j);
        uint32_t span = new_max - new_min + 1;

        if (span <= mSpanBlocks) {
          if (L.used_bins.empty()) {
            L.used_bins.assign(mNumSlots, 0);
          }
          L.used_bins[r] = 1;
          L.min_block = new_min;
          L.max_block = new_max;
          mItemToLayerIdx[idx] = li;
          placed = true;
          break;
        }
      }

      if (!placed) {
        Layer nl;
        nl.min_block = j;
        nl.max_block = j;
        nl.used_bins.assign(mNumSlots, 0);
        nl.used_bins[r] = 1;
        layers.push_back(std::move(nl));
        mItemToLayerIdx[idx] = layers.size() - 1;
      }
    }
  }

  mNumLayers = layers.size();
  mLayerMinBlock.resize(mNumLayers);
  mLayerMaxBlock.resize(mNumLayers);
  for (uint32_t l = 0; l < mNumLayers; ++l) {
    mLayerMinBlock[l] = layers[l].min_block;
    mLayerMaxBlock[l] = layers[l].max_block;
  }

  std::vector<std::vector<uint32_t>> bin_layers(mNumSlots);
  mLayerBins.resize(mNumLayers);
  for (uint32_t l = 0; l < mNumLayers; l++) {
    mLayerBins[l].assign(mNumSlots, UINT32_MAX);
  }

  for (uint32_t i = 0; i < mN; ++i) {
    uint32_t bin = item_binidx[i];
    uint32_t l = mItemToLayerIdx[i];
    bin_layers[bin].push_back(l);
    mLayerBins[l][bin] = i;
  }

  bin_sizes.resize(mNumSlots);

  for (uint32_t k = 0; k < mNumSlots; ++k) {
    auto &nonempty_layers = bin_layers[k];
    if (nonempty_layers.empty()) {
      continue;
    }
    // last (nonempty) layer of k-th bin
    uint32_t last_layer = 0;
    for (auto &layeridx : nonempty_layers) {
      last_layer = max(layeridx, last_layer);
    }
    bin_sizes[k] = last_layer + 1;
    BitVector oc_indicator(last_layer + 1);
    for (auto &layeridx : nonempty_layers) {
      oc_indicator[layeridx] = 1;
    }
    occupy_indicator_flat.append(oc_indicator);
  }

  ot_idx.resize(mN);
  size_t idx = 0;

  if (mRpmt) {
    for (size_t lay = 0; lay < mNumLayers; ++lay) {
      for (uint32_t bin = 0; bin < mNumSlots; ++bin) {
        if (mLayerBins[lay][bin] != UINT32_MAX) {
          ot_idx[idx++] = mLayerBins[lay][bin];
        }
      }
    }
  } else {
    maskings.resize(mN);
    ptxts_mask.resize(mNumLayers);
    vector<vector<uint64_t>> raw_masks(mNumLayers);
    for (size_t lay = 0; lay < mNumLayers; lay++) {
      raw_masks[lay].resize(mNumSlots);
      mPrng.get<uint64_t>(raw_masks[lay]);
      for (uint32_t bin = 0; bin < mNumSlots; bin++) {
        raw_masks[lay][bin] =
            seal::util::barrett_reduce_64(raw_masks[lay][bin], mModulus);
        if (mLayerBins[lay][bin] != UINT32_MAX) {
          ot_idx[idx] = mLayerBins[lay][bin];
          maskings[idx++] = raw_masks[lay][bin];
        }
      }
      mBatchEncoder->encode(raw_masks[lay], ptxts_mask[lay]);
    }
  }
}

void SspmtSender::init(uint32_t n, uint32_t nReceiver, sspmtParams ssParams,
                       oc::block seed) {
  EncryptionParameters parms(scheme_type::bgv);
  mNumSlots = ssParams.heNumSlots;
  parms.set_poly_modulus_degree(mNumSlots);

  mN = n;
  mNreceiver = nReceiver;
  mM = roundUpTo(ssParams.bandExpansion * n, mNumSlots);
  mW = ssParams.bandWidth;
  mNumBatch = mM / mNumSlots;
  mWrap = divCeil(mW * mNumSlots, mM) + 1;
  mSpanBlocks = ssParams.span_blocks;
  mPrng.SetSeed(seed);

  mItemToBlockIdx.resize(mN);
  mItemToLayerIdx.resize(mN);

  parms.set_coeff_modulus(
      CoeffModulus::Create(mNumSlots, ssParams.heCoeffModulus));
  mModulus = PlainModulus::Batching(mNumSlots, ssParams.hePlainModulusBits);
  parms.set_plain_modulus(mModulus);

  mContext = make_shared<SEALContext>(parms);
  mBatchEncoder = make_unique<BatchEncoder>(*mContext);
  mEvaluator = make_unique<Evaluator>(*mContext);
};

Proto SspmtSender::run(const std::vector<oc::block> &Y, Socket &chl) {
  preprocess(Y);

  SEALContext context = *mContext;

  string recvstring;
  co_await chl.recvResize(recvstring);
  cout << "Sender receives " << mNumBatch << " X " << mWrap
       << " Ctxts of OKVS Encoding, " << recvstring.size() << " Bytes" << endl;

  stringstream recvstream;
  recvstream.write(recvstring.data(), recvstring.size());

  vector<vector<Ciphertext>> encoded_in_he(mWrap);
  for (size_t i = 0; i < mWrap; i++) {
    encoded_in_he[i].resize(mNumBatch);
    for (size_t j = 0; j < mNumBatch; j++) {
      encoded_in_he[i][j].unsafe_load(context, recvstream);
    }
  }

  setTimePoint("Sender::Recv ctxts & Serialize");

  vector<Ciphertext> decoded_in_he;

  encrypted_decode(encoded_in_he, decoded_in_he);

  stringstream sendstream;
  for (size_t i = 0; i < decoded_in_he.size(); i++) {
    decoded_in_he[i].save(sendstream);
  }

  cout << "Sender sends " << decoded_in_he.size() << " decoded ctxts, "
       << sendstream.str().size() << " Bytes" << endl;

  co_await chl.send(decoded_in_he.size());
  co_await chl.send(move(sendstream.str()));
  co_await chl.send(move(bin_sizes));
  co_await chl.send(move(occupy_indicator_flat));

  setTimePoint("Sender::Serialize & Send");
}

Proto SspmtSender::run(const std::vector<oc::block> &Y, oc::BitVector &results,
                       Socket &chl) {
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

void SspmtSender::preprocess(const std::vector<oc::block> &Y) {
  PrimeFieldOkvs okvs;
  // okvs.setTimer(getTimer());
  okvs.init(Y.size(), mM, mW, mModulus);
  vector<uint64_t> bands_flat(mN * mW);
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

  sequencing(start_pos_spacing);

  setTimePoint("Sender::Sequencing");

  const uint64_t *bands = bands_flat.data();

  ptxts_diags.resize(mWrap);
  for (size_t k = 0; k < mWrap; ++k) {
    ptxts_diags[k].assign(mNumLayers * mNumBatch, seal::Plaintext(mNumSlots));
  }

  std::vector<uint64_t> plainVec(mNumSlots);

  for (uint32_t i = 0; i < mNumLayers; ++i) {
    uint32_t Bmin = mLayerMinBlock[i];
    uint32_t Bmax = mLayerMaxBlock[i] + (mW - 1);

    for (uint32_t B = Bmin; B <= Bmax; ++B) {
      uint32_t k = B / mNumBatch;  // diag index
      uint32_t j = B % mNumBatch;  // batch (block) index
      if (k >= mWrap) continue;

      size_t outIdx = i * mNumBatch + j;

      bool all_zero = true;
      for (uint32_t bin = 0; bin < mNumSlots; ++bin) {
        uint32_t item = mLayerBins[i][bin];
        uint64_t val = 0;
        if (item != UINT32_MAX) {
          uint32_t startBlk = mItemToBlockIdx[item];
          if (B >= startBlk && B < startBlk + mW) {
            uint32_t w = B - startBlk;
            val = bands[item * mW + w];
          }
        }
        plainVec[bin] = val;
        if (val != 0) all_zero = false;
      }

      if (!all_zero) {
        mBatchEncoder->encode(plainVec, ptxts_diags[k][outIdx]);
      }
    }
  }
  setTimePoint("Sender::HE Encode");
}

void SspmtSender::encrypted_decode(
    const std::vector<std::vector<seal::Ciphertext>> &encoded_in_he,
    std::vector<seal::Ciphertext> &decoded_in_he) {
  decoded_in_he.resize(mNumLayers);

  for (size_t i = 0; i < mNumLayers; ++i) {
    bool initialized = false;
    Ciphertext tmp;
    for (size_t j = 0; j < mNumBatch; ++j) {
      const size_t idx = i * mNumBatch + j;
      for (size_t k = 0; k < mWrap; k++) {
        const Plaintext &pt_diag = ptxts_diags[k][idx];
        if (!pt_diag.is_zero()) {
          if (!initialized) {
            mEvaluator->multiply_plain(encoded_in_he[k][j], pt_diag,
                                       decoded_in_he[i]);
            initialized = true;
          } else {
            mEvaluator->multiply_plain(encoded_in_he[k][j], pt_diag, tmp);
            mEvaluator->add_inplace(decoded_in_he[i], tmp);
          }
        }
      }
    }
  }

  if (!mRpmt) {
    for (size_t i = 0; i < mNumLayers; ++i) {
      mEvaluator->add_plain_inplace(decoded_in_he[i], ptxts_mask[i]);
    }
  }

  for (size_t i = 0; i < decoded_in_he.size(); i++) {
    mEvaluator->mod_switch_to_next_inplace(decoded_in_he[i]);
  }
  setTimePoint("Sender::Encrypted OKVS Decoding");
}

void SspmtReceiver::init(uint32_t n, uint32_t nSender, sspmtParams ssParams,
                         oc::block seed) {
  EncryptionParameters parms(scheme_type::bgv);
  mNumSlots = ssParams.heNumSlots;
  parms.set_poly_modulus_degree(mNumSlots);

  mN = n;
  mNsender = nSender;
  mM = roundUpTo(ssParams.bandExpansion * n, mNumSlots);
  mW = ssParams.bandWidth;
  mNumBatch = mM / mNumSlots;
  mWrap = divCeil(mW * mNumSlots, mM) + 1;
  mPrng.SetSeed(seed);

  parms.set_coeff_modulus(
      CoeffModulus::Create(mNumSlots, ssParams.heCoeffModulus));
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

Proto SspmtReceiver::run(const std::vector<oc::block> &X,
                         oc::BitVector &results, Socket &chl) {
  stringstream sendstream;
  encode_and_encrypt(X, sendstream);
  co_await chl.send(move(sendstream.str()));
  setTimePoint("Receiver::Serialize & Send");

  size_t decoded_he_size;
  string recvstring;
  co_await chl.recv(decoded_he_size);
  co_await chl.recvResize(recvstring);
  co_await chl.recvResize(bin_sizes);
  auto len = 0;
  for (size_t i = 0; i < mNumSlots; i++) {
    len += bin_sizes[i];
  }
  oc::BitVector occupy_indicator_flat(len);
  co_await chl.recv(occupy_indicator_flat);
  auto offset = 0;
  occupy_indicator.resize(mNumSlots);
  for (size_t i = 0; i < mNumSlots; i++) {
    occupy_indicator[i].append(occupy_indicator_flat.data(), bin_sizes[i],
                               offset);
    offset += bin_sizes[i];
  }

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
    for (size_t i = 0; i < mNsender; i++) {
      results[i] = (dec_results[i] == mIndicatorStr) ? 1 : 0;
    }
  } else {
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

void SspmtReceiver::encode_and_encrypt(const std::vector<oc::block> &X,
                                       stringstream &ctxtstream) {
  mIndicatorStr =
      seal::util::barrett_reduce_64(mPrng.get<uint64_t>(), mModulus);
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

  for (size_t k = 0; k < mWrap; k++) {
    auto lastIdx = (k > 0) ? mNumBatch - 1 : mNumBatch;
    for (uint32_t i = 0; i < lastIdx; i++) {
      std::copy_n(&encoded[i * mNumSlots + k], mNumSlots, plainVec.data());
      mBatchEncoder->encode(plainVec, ptxt);
      mEncryptor->encrypt_symmetric(ptxt).save(ctxtstream);
    }
    if (k > 0) {
      std::copy_n(&encoded[lastIdx * mNumSlots + k], mNumSlots - k,
                  plainVec.data());
      std::copy_n(&encoded[0], k, plainVec.data() + mNumSlots - k);
      mBatchEncoder->encode(plainVec, ptxt);
      mEncryptor->encrypt_symmetric(ptxt).save(ctxtstream);
    }
  }
  setTimePoint("Receiver::Encryption");
}

void SspmtReceiver::decrypt(const std::vector<seal::Ciphertext> &decoded_in_he,
                            std::vector<uint64_t> &dec_results) {
  const size_t L = decoded_in_he.size();
  vector<Plaintext> ptxts(L);

  //   cout << "Noise Budget: "
  //    << mDecryptor->invariant_noise_budget(decoded_in_he[0]) << endl;

  for (size_t i = 0; i < L; i++) {
    mDecryptor->decrypt(decoded_in_he[i], ptxts[i]);
  }

  dec_results.resize(mNsender);
  vector<uint64_t> decodeVec(mNumSlots);

  auto idx = 0;
  for (size_t i = 0; i < L; i++) {
    mBatchEncoder->decode(ptxts[i], decodeVec);

    for (size_t j = 0; j < mNumSlots; j++) {
      if (i < occupy_indicator[j].size()) {
        if (occupy_indicator[j][i]) {
          dec_results[idx++] = decodeVec[j];
        }
      }
    }
  }
  setTimePoint("Receiver::Decrypt");
}
}  // namespace rlweOkvs