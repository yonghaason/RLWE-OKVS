#include "rpmt.h"

#include <cstdlib>
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

namespace rlweOkvs
{
  namespace
  {
    uint32_t chunkEnvOrDefault(const char *name, uint32_t fallback)
    {
      if (const char *value = std::getenv(name))
      {
        auto parsed = static_cast<uint32_t>(std::strtoul(value, nullptr, 10));
        if (parsed > 0)
        {
          return parsed;
        }
      }
      return fallback;
    }

    uint32_t encodedBatchChunk()
    {
      return chunkEnvOrDefault("RLWE_OKVS_ENCODED_BATCH_CHUNK", 8);
    }

    uint32_t decodedLayerChunk()
    {
      return chunkEnvOrDefault("RLWE_OKVS_DECODED_LAYER_CHUNK", 8);
    }

    uint64_t totalEncodedCipherCount(uint32_t numBatch, uint32_t width, uint32_t maxWrap)
    {
      (void)width;
      (void)maxWrap;
      return numBatch;
    }

    uint32_t halfWrapCount(uint32_t width, uint32_t numHalfBatch)
    {
      if (width == 0)
      {
        return 1;
      }
      return 1 + divCeil(width, numHalfBatch);
    }

    uint32_t activeShiftBatchCount(
        uint32_t carry, uint32_t width, uint32_t numHalfBatch, uint32_t numBatch)
    {
      if (carry == 0)
      {
        return 0;
      }

      const int64_t activeHalfBlocks =
          static_cast<int64_t>(width) - 1 - static_cast<int64_t>(carry - 1) * numHalfBatch;
      if (activeHalfBlocks <= 0)
      {
        return 0;
      }

      return std::min<uint32_t>(
          numBatch, static_cast<uint32_t>(divCeil(static_cast<uint32_t>(activeHalfBlocks), 2u)));
    }

  }

  void RpmtSender::sequencing(const std::vector<uint32_t> &start_pos_spacing)
  {
    occupy_indicator_flat = oc::BitVector();
    std::vector<uint32_t> item_binidx(mN);
    uint32_t max_block = 0;

    for (uint32_t i = 0; i < mN; ++i)
    {
      uint32_t pos = start_pos_spacing[i];
      uint32_t bin = pos % mHalfSlots;
      uint32_t blk = pos / mHalfSlots;

      item_binidx[i] = bin;
      mItemToBlockIdx[i] = blk;
      if (blk > max_block)
        max_block = blk;
    }

    std::vector<std::vector<uint32_t>> block_items(max_block + 1);
    for (uint32_t i = 0; i < mN; ++i)
    {
      block_items[mItemToBlockIdx[i]].push_back(i);
    }

    struct Layer
    {
      uint32_t min_block;
      uint32_t max_block;
      std::vector<uint8_t> used_slots;
    };
    std::vector<Layer> layers;

    for (uint32_t blk = 0; blk <= max_block; ++blk)
    {
      auto &bucket = block_items[blk];
      for (uint32_t idx : bucket)
      {
        uint32_t r = item_binidx[idx];     // half-row residue
        uint32_t j = mItemToBlockIdx[idx]; // half-block
        bool placed = false;
        for (uint32_t li = 0; li < layers.size(); ++li)
        {
          Layer &L = layers[li];
          const uint32_t top_slot = r;
          const uint32_t bottom_slot = mHalfSlots + r;

          if (L.used_slots[top_slot] && L.used_slots[bottom_slot])
          {
            continue;
          }

          uint32_t new_min = std::min(L.min_block, j);
          uint32_t new_max = std::max(L.max_block, j);
          uint32_t span = new_max - new_min + 1;

          if (span <= mSpanBlocks)
          {
            const uint32_t slot = L.used_slots[top_slot] ? bottom_slot : top_slot;
            L.used_slots[slot] = 1;
            L.min_block = new_min;
            L.max_block = new_max;
            mItemToLayerIdx[idx] = li;
            mItemToSlotIdx[idx] = slot;
            placed = true;
            break;
          }
        }

        if (!placed)
        {
          Layer nl;
          nl.min_block = j;
          nl.max_block = j;
          nl.used_slots.assign(mNumSlots, 0);
          nl.used_slots[r] = 1;
          layers.push_back(std::move(nl));
          mItemToLayerIdx[idx] = layers.size() - 1;
          mItemToSlotIdx[idx] = r;
        }
      }
    }

    mNumLayers = layers.size();
    mLayerMinBlock.resize(mNumLayers);
    mLayerMaxBlock.resize(mNumLayers);
    for (uint32_t l = 0; l < mNumLayers; ++l)
    {
      mLayerMinBlock[l] = layers[l].min_block;
      mLayerMaxBlock[l] = layers[l].max_block;
    }

    std::vector<std::vector<uint32_t>> slot_layers(mNumSlots);
    mLayerBins.resize(mNumLayers);
    for (uint32_t l = 0; l < mNumLayers; l++)
    {
      mLayerBins[l].assign(mNumSlots, UINT32_MAX);
    }

    for (uint32_t i = 0; i < mN; ++i)
    {
      uint32_t slot = mItemToSlotIdx[i];
      uint32_t l = mItemToLayerIdx[i];
      slot_layers[slot].push_back(l);
      mLayerBins[l][slot] = i;
    }

    last_layer_per_bin.resize(mNumSlots);

    for (uint32_t k = 0; k < mNumSlots; ++k)
    {
      auto &nonempty_layers = slot_layers[k];
      if (nonempty_layers.empty())
      {
        last_layer_per_bin[k] = 0;
        continue;
      }

      uint32_t last_layer = 0;
      for (auto &layeridx : nonempty_layers)
      {
        last_layer = max(layeridx, last_layer);
      }
      last_layer_per_bin[k] = last_layer + 1;
      BitVector oc_indicator(last_layer + 1);
      for (auto &layeridx : nonempty_layers)
      {
        oc_indicator[layeridx] = 1;
      }
      occupy_indicator_flat.append(oc_indicator);
    }

    ot_idx.resize(mN);
    size_t idx = 0;

    if (!mSharedOutput)
    {
      for (size_t lay = 0; lay < mNumLayers; ++lay)
      {
        for (uint32_t bin = 0; bin < mNumSlots; ++bin)
        {
          if (mLayerBins[lay][bin] != UINT32_MAX)
          {
            ot_idx[idx++] = mLayerBins[lay][bin];
          }
        }
      }
    }
    else
    {
      maskings.resize(mN);
      ptxts_mask.resize(mNumLayers);
      vector<vector<uint64_t>> raw_masks(mNumLayers);
      for (size_t lay = 0; lay < mNumLayers; lay++)
      {
        raw_masks[lay].resize(mNumSlots);
        mPrng.get<uint64_t>(raw_masks[lay]);
        for (uint32_t bin = 0; bin < mNumSlots; bin++)
        {
          raw_masks[lay][bin] =
              seal::util::barrett_reduce_64(raw_masks[lay][bin], mModulus);
          if (mLayerBins[lay][bin] != UINT32_MAX)
          {
            ot_idx[idx] = mLayerBins[lay][bin];
            maskings[idx++] = raw_masks[lay][bin];
          }
        }
        mBatchEncoder->encode(raw_masks[lay], ptxts_mask[lay]);
      }
    }
  }

  void RpmtSender::init(uint32_t n, uint32_t nReceiver, rpmtParams ssParams,
                         oc::block seed)
  {
    EncryptionParameters parms(scheme_type::bgv);
    mNumSlots = ssParams.heNumSlots;
    parms.set_poly_modulus_degree(mNumSlots);

    mN = n;
    mNreceiver = nReceiver;
    mM = roundUpTo(ssParams.bandExpansion * n, mNumSlots);
    mW = ssParams.bandWidth;
    mHalfSlots = static_cast<uint32_t>(mNumSlots / 2);
    mNumBatch = mM / mNumSlots;
    mNumHalfBatch = mM / mHalfSlots;
    mWrap = halfWrapCount(mW, mNumHalfBatch);
    mSpanBlocks = ssParams.span_blocks;
    mPrng.SetSeed(seed);

    mItemToBlockIdx.resize(mN);
    mItemToLayerIdx.resize(mN);
    mItemToSlotIdx.resize(mN);

    parms.set_coeff_modulus(
        CoeffModulus::Create(mNumSlots, ssParams.heCoeffModulus));
    mModulus = seal::util::get_primes(2*mNumSlots, ssParams.hePlainModulusBits, 4)[3];
    // mModulus = PlainModulus::Batching(mNumSlots, ssParams.hePlainModulusBits);
    parms.set_plain_modulus(mModulus);

    mContext = make_shared<SEALContext>(parms, true, sec_level_type::none);
    mBatchEncoder = make_unique<BatchEncoder>(*mContext);
    mEvaluator = make_unique<Evaluator>(*mContext);
  };

  Proto RpmtSender::run(const std::vector<oc::block> &Y, Socket &chl)
  {
    std::vector<oc::block> FY;
    mOprfSender.init(mN, mNreceiver, mPrng.get());
    mOprfSender.setTimer(getTimer());
    co_await mOprfSender.run(Y, FY, chl);
    setTimePoint("Sender::OPRF");

    preprocess(FY);

    string gk_string;
    co_await chl.recvResize(gk_string);
    cout << "Sender receives Galois keys, " << gk_string.size() << " Bytes" << endl;
    stringstream gk_stream;
    gk_stream.str(gk_string);
    mGaloisKeys.load(*mContext, gk_stream);
    setTimePoint("Sender::Recv Galois Keys");

    vector<Ciphertext> encoded_in_he(mNumBatch);
    co_await recv_encoded_chunks(encoded_in_he, chl);
    co_await send_decoded_chunks(encoded_in_he, chl);
  }

  Proto RpmtSender::run(const std::vector<oc::block> &Y, oc::BitVector &results,
                         Socket &chl)
  {
    assert(mSharedOutput && "Sender obtains result only when ssPMT.");
    co_await run(Y, chl);

    // GMW with maskings
    u64 keyBitLength = 40 + oc::log2ceil(Y.size());
    u64 keyByteLength = oc::divCeil(keyBitLength, 8);

    oc::Matrix<u8> gmwin;
    gmwin.resize(Y.size(), keyByteLength, oc::AllocType::Uninitialized);
    for (size_t i = 0; i < Y.size(); i++)
    {
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

  void RpmtSender::preprocess(const std::vector<oc::block> &Y)
  {
    PrimeFieldOkvs okvs;
    // okvs.setTimer(getTimer());

    okvs.init(Y.size(), mM, mW, mModulus);
    vector<uint64_t> bands_flat(mN * mW);
    vector<uint32_t> start_pos(mN);
    okvs.generate_band(Y, bands_flat, start_pos, oc::ZeroBlock);

    setTimePoint("Sender::Generate Bands");

    vector<uint32_t> start_pos_spacing(mN);
    for (uint32_t i = 0; i < mN; i++)
    {
      auto position = start_pos[i];
      uint32_t q = position / mNumHalfBatch;
      uint32_t r = position % mNumHalfBatch;
      start_pos_spacing[i] = r * mHalfSlots + q;
    }

    sequencing(start_pos_spacing);

    setTimePoint("Sender::Sequencing");

    struct SlotMeta
    {
      uint32_t slot;
      const uint64_t *band_ptr;
      uint32_t start_blk;
    };

    const uint64_t *bands = bands_flat.data();

    ptxts_diags.clear();
    ptxts_diags_swapped.clear();
    ptxts_diags_shifted.clear();
    ptxts_diags_shifted_swapped.clear();
    ptxts_diags.resize(mNumLayers);
    ptxts_diags_swapped.resize(mNumLayers);
    ptxts_diags_shifted.resize(mWrap > 1 ? mWrap - 1 : 0);
    ptxts_diags_shifted_swapped.resize(mWrap > 1 ? mWrap - 1 : 0);
    for (size_t i = 0; i < mNumLayers; ++i)
    {
      ptxts_diags[i].resize(mNumBatch);
      ptxts_diags_swapped[i].resize(mNumBatch);
    }
    for (uint32_t carry = 1; carry < mWrap; ++carry)
    {
      ptxts_diags_shifted[carry - 1].resize(mNumLayers);
      ptxts_diags_shifted_swapped[carry - 1].resize(mNumLayers);
      for (size_t i = 0; i < mNumLayers; ++i)
      {
        ptxts_diags_shifted[carry - 1][i].resize(mNumBatch);
        ptxts_diags_shifted_swapped[carry - 1][i].resize(mNumBatch);
      }
    }

    std::vector<SlotMeta> layer_meta;
    layer_meta.reserve(mNumSlots);
    std::vector<uint64_t> mainVec(mNumSlots, 0);
    std::vector<uint64_t> swappedVec(mNumSlots, 0);
    using Entry = std::pair<uint32_t, uint64_t>;
    std::vector<std::vector<Entry>> mainEntries(mNumBatch);
    std::vector<std::vector<Entry>> swappedEntries(mNumBatch);
    std::vector<std::vector<std::vector<Entry>>> shiftedEntries(
        mWrap > 1 ? mWrap - 1 : 0, std::vector<std::vector<Entry>>(mNumBatch));
    std::vector<std::vector<std::vector<Entry>>> shiftedSwappedEntries(
        mWrap > 1 ? mWrap - 1 : 0, std::vector<std::vector<Entry>>(mNumBatch));

    for (uint32_t i = 0; i < mNumLayers; ++i)
    {
      layer_meta.clear();
      for (uint32_t batch = 0; batch < mNumBatch; ++batch)
      {
        mainEntries[batch].clear();
        swappedEntries[batch].clear();
      }
      for (uint32_t carry = 1; carry < mWrap; ++carry)
      {
        const uint32_t activeBatchCount =
            activeShiftBatchCount(carry, mW, mNumHalfBatch, mNumBatch);
        for (uint32_t batch = 0; batch < activeBatchCount; ++batch)
        {
          shiftedEntries[carry - 1][batch].clear();
          shiftedSwappedEntries[carry - 1][batch].clear();
        }
      }

      for (uint32_t slot = 0; slot < mNumSlots; ++slot)
      {
        uint32_t item = mLayerBins[i][slot];
        if (item == UINT32_MAX)
        {
          continue;
        }

        layer_meta.push_back(SlotMeta{
            slot,
            bands + static_cast<size_t>(item) * mW,
            mItemToBlockIdx[item]});
      }

      for (const auto &sm : layer_meta)
      {
        const bool isTop = sm.slot < mHalfSlots;
        const uint32_t rowBase = isTop ? 0 : mHalfSlots;
        const uint32_t logicalBin = sm.slot - rowBase;
        for (uint32_t w = 0; w < mW; ++w)
        {
          const uint32_t totalHalfBlock = sm.start_blk + w;
          const uint32_t blockCarry = totalHalfBlock / mNumHalfBatch;
          const uint32_t halfBlock = totalHalfBlock - blockCarry * mNumHalfBatch;
          const uint32_t batch = halfBlock / 2;
          const uint64_t coeff = sm.band_ptr[w];
          auto &entryList =
              blockCarry == 0
                  ? (((halfBlock & 1U) == 0)
                         ? (isTop ? mainEntries[batch] : swappedEntries[batch])
                         : (isTop ? swappedEntries[batch] : mainEntries[batch]))
                  : (((halfBlock & 1U) == 0)
                         ? (isTop ? shiftedEntries[blockCarry - 1][batch]
                                  : shiftedSwappedEntries[blockCarry - 1][batch])
                         : (isTop ? shiftedSwappedEntries[blockCarry - 1][batch]
                                  : shiftedEntries[blockCarry - 1][batch]));

          entryList.push_back({rowBase + logicalBin, coeff});
        }
      }

      for (uint32_t batch = 0; batch < mNumBatch; ++batch)
      {
        if (!mainEntries[batch].empty())
        {
          for (const auto &[slot, coeff] : mainEntries[batch])
          {
            mainVec[slot] = coeff;
          }
          mBatchEncoder->encode(mainVec, ptxts_diags[i][batch]);
          for (const auto &[slot, coeff] : mainEntries[batch])
          {
            (void)coeff;
            mainVec[slot] = 0;
          }
        }
        if (!swappedEntries[batch].empty())
        {
          for (const auto &[slot, coeff] : swappedEntries[batch])
          {
            swappedVec[slot] = coeff;
          }
          mBatchEncoder->encode(swappedVec, ptxts_diags_swapped[i][batch]);
          for (const auto &[slot, coeff] : swappedEntries[batch])
          {
            (void)coeff;
            swappedVec[slot] = 0;
          }
        }
      }

      for (uint32_t carry = 1; carry < mWrap; ++carry)
      {
        const uint32_t activeBatchCount =
            activeShiftBatchCount(carry, mW, mNumHalfBatch, mNumBatch);
        for (uint32_t batch = 0; batch < activeBatchCount; ++batch)
        {
          if (!shiftedEntries[carry - 1][batch].empty())
          {
            for (const auto &[slot, coeff] : shiftedEntries[carry - 1][batch])
            {
              mainVec[slot] = coeff;
            }
            mBatchEncoder->encode(mainVec, ptxts_diags_shifted[carry - 1][i][batch]);
            for (const auto &[slot, coeff] : shiftedEntries[carry - 1][batch])
            {
              (void)coeff;
              mainVec[slot] = 0;
            }
          }
          if (!shiftedSwappedEntries[carry - 1][batch].empty())
          {
            for (const auto &[slot, coeff] : shiftedSwappedEntries[carry - 1][batch])
            {
              swappedVec[slot] = coeff;
            }
            mBatchEncoder->encode(
                swappedVec, ptxts_diags_shifted_swapped[carry - 1][i][batch]);
            for (const auto &[slot, coeff] : shiftedSwappedEntries[carry - 1][batch])
            {
              (void)coeff;
              swappedVec[slot] = 0;
            }
          }
        }
      }
    }

    setTimePoint("Sender::HE Encode");
  }

  Proto RpmtSender::recv_encoded_chunks(
      std::vector<seal::Ciphertext> &encoded_in_he,
      Socket &chl)
  {
    SEALContext context = *mContext;    
    size_t recvBytes = 0;
    string recvstring;
    stringstream recvstream;

    const uint32_t chunkSize = encodedBatchChunk();
    for (uint32_t jBegin = 0; jBegin < mNumBatch; jBegin += chunkSize)
    {
      const uint32_t jEnd = std::min<uint32_t>(mNumBatch, jBegin + chunkSize);
      co_await chl.recvResize(recvstring);
      recvBytes += recvstring.size();

      recvstream.clear();
      recvstream.str(recvstring);

      for (uint32_t j = jBegin; j < jEnd; ++j)
      {
        encoded_in_he[j].unsafe_load(context, recvstream);
      }
    }

    cout << "Sender receives " << mNumBatch
         << " Ctxts of OKVS Encoding, " << recvBytes << " Bytes" << endl;

    setTimePoint("Sender::Recv ctxts & Serialize");
  }

  Proto RpmtSender::send_decoded_chunks(
      const std::vector<seal::Ciphertext> &encoded_in_he,
      Socket &chl)
  {
    co_await chl.send(mNumLayers);

    size_t sentBytes = 0;
    std::vector<Ciphertext> decoded_chunk;
    std::vector<Ciphertext> swapped_rows;
    std::vector<std::vector<Ciphertext>> shifted_rows(mWrap > 1 ? mWrap - 1 : 0);
    std::vector<std::vector<Ciphertext>> shifted_swapped_rows(mWrap > 1 ? mWrap - 1 : 0);
    stringstream sendstream;

    swapped_rows.resize(mNumBatch);
    for (uint32_t j = 0; j < mNumBatch; ++j)
    {
      mEvaluator->rotate_columns(encoded_in_he[j], mGaloisKeys, swapped_rows[j]);
    }
    for (uint32_t carry = 1; carry < mWrap; ++carry)
    {
      const uint32_t activeBatchCount =
          activeShiftBatchCount(carry, mW, mNumHalfBatch, mNumBatch);
      shifted_rows[carry - 1].resize(activeBatchCount);
      shifted_swapped_rows[carry - 1].resize(activeBatchCount);
      for (uint32_t j = 0; j < activeBatchCount; ++j)
      {
        mEvaluator->rotate_rows(
            encoded_in_he[j], static_cast<int>(carry), mGaloisKeys,
            shifted_rows[carry - 1][j]);
        mEvaluator->rotate_columns(
            shifted_rows[carry - 1][j], mGaloisKeys, shifted_swapped_rows[carry - 1][j]);
      }
    }
    setTimePoint("Sender::Precompute Rotations");

    const uint32_t chunkSize = decodedLayerChunk();
    for (uint32_t layerBegin = 0; layerBegin < mNumLayers; layerBegin += chunkSize)
    {
      const uint32_t layerEnd =
          std::min<uint32_t>(mNumLayers, layerBegin + chunkSize);

      encrypted_decode(encoded_in_he, swapped_rows, shifted_rows,
                       shifted_swapped_rows, decoded_chunk,
                       layerBegin, layerEnd);

      sendstream.clear();
      sendstream.str("");
      for (size_t i = 0; i < decoded_chunk.size(); ++i)
      {
        decoded_chunk[i].save(sendstream);
      }

      auto payload = sendstream.str();
      sentBytes += payload.size();
      co_await chl.send(move(payload));
    }

    cout << "Sender sends " << mNumLayers << " decoded ctxts, "
         << sentBytes << " Bytes" << endl;

    co_await chl.send(move(last_layer_per_bin));
    co_await chl.send(move(occupy_indicator_flat));

    setTimePoint("Sender::Encrypted OKVS Decoding & Send Back");
  }

  void RpmtSender::encrypted_decode(
      const std::vector<seal::Ciphertext> &encoded_in_he,
      const std::vector<seal::Ciphertext> &swapped_rows,
      const std::vector<std::vector<seal::Ciphertext>> &shifted_rows,
      const std::vector<std::vector<seal::Ciphertext>> &shifted_swapped_rows,
      std::vector<seal::Ciphertext> &decoded_in_he,
      uint32_t layerBegin,
      uint32_t layerEnd)
  {
    decoded_in_he.resize(layerEnd - layerBegin);

    for (uint32_t i = layerBegin; i < layerEnd; ++i)
    {
      bool initialized = false;
      Ciphertext &out = decoded_in_he[i - layerBegin];

      for (uint32_t j = 0; j < mNumBatch; ++j)
      {
        Ciphertext prod;

        auto accumulate = [&](const Ciphertext &src, const Plaintext &plain)
        {
          if (plain.is_zero())
          {
            return;
          }

          mEvaluator->multiply_plain(src, plain, prod);
          if (!initialized)
          {
            out = prod;
            initialized = true;
          }
          else
          {
            mEvaluator->add_inplace(out, prod);
          }
        };

        accumulate(encoded_in_he[j], ptxts_diags[i][j]);
        accumulate(swapped_rows[j], ptxts_diags_swapped[i][j]);
        for (uint32_t carry = 1; carry < mWrap; ++carry)
        {
          if (j < shifted_rows[carry - 1].size() &&
              !ptxts_diags_shifted[carry - 1][i][j].is_zero())
          {
            accumulate(shifted_rows[carry - 1][j], ptxts_diags_shifted[carry - 1][i][j]);
          }
          if (j < shifted_swapped_rows[carry - 1].size() &&
              !ptxts_diags_shifted_swapped[carry - 1][i][j].is_zero())
          {
            accumulate(
                shifted_swapped_rows[carry - 1][j],
                ptxts_diags_shifted_swapped[carry - 1][i][j]);
          }
        }
      }

      if (!initialized)
      {
        out = encoded_in_he[0];
        mEvaluator->sub_inplace(out, encoded_in_he[0]);
      }

      if (mSharedOutput)
      {
        mEvaluator->add_plain_inplace(out, ptxts_mask[i]);
      }

      mEvaluator->mod_switch_to_next_inplace(out);
    }
  }

  void RpmtReceiver::init(uint32_t n, uint32_t nSender, rpmtParams ssParams,
                           oc::block seed)
  {
    EncryptionParameters parms(scheme_type::bgv);
    mNumSlots = ssParams.heNumSlots;
    parms.set_poly_modulus_degree(mNumSlots);

    mN = n;
    mNsender = nSender;
    mM = roundUpTo(ssParams.bandExpansion * n, mNumSlots);
    mW = ssParams.bandWidth;
    mHalfSlots = static_cast<uint32_t>(mNumSlots / 2);
    mNumBatch = mM / mNumSlots;
    mNumHalfBatch = mM / mHalfSlots;
    mWrap = halfWrapCount(mW, mNumHalfBatch);
    mPrng.SetSeed(seed);

    parms.set_coeff_modulus(
        CoeffModulus::Create(mNumSlots, ssParams.heCoeffModulus));

    mModulus = seal::util::get_primes(2*mNumSlots, ssParams.hePlainModulusBits, 4)[3];
    // mModulus = PlainModulus::Batching(mNumSlots, ssParams.hePlainModulusBits);
    parms.set_plain_modulus(mModulus);

    mContext = make_shared<SEALContext>(parms, true, sec_level_type::none);
    KeyGenerator keygen(*mContext);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    std::vector<int> galois_steps;
    galois_steps.push_back(0);
    for (uint32_t step = 1; step < mWrap; ++step)
    {
      galois_steps.push_back(static_cast<int>(step));
    }
    keygen.create_galois_keys(galois_steps, mGaloisKeys);

    mEncryptor = make_unique<Encryptor>(*mContext, public_key);
    mEncryptor->set_secret_key(secret_key);

    mBatchEncoder = make_unique<BatchEncoder>(*mContext);
    mDecryptor = make_unique<Decryptor>(*mContext, secret_key);
  };

  Proto RpmtReceiver::run(const std::vector<oc::block> &X,
                           oc::BitVector &results, Socket &chl)
  {
    std::vector<oc::block> FX;
    mOprfReceiver.init(mN, mNsender, mPrng.get());
    mOprfReceiver.setTimer(getTimer());
    co_await mOprfReceiver.run(X, FX, chl);
    setTimePoint("Receiver::OPRF");

    stringstream gk_stream;
    mGaloisKeys.save(gk_stream);
    auto gk_payload = gk_stream.str();
    cout << "Receiver sends Galois keys, " << gk_payload.size() << " Bytes" << endl;
    co_await chl.send(std::move(gk_payload));
    setTimePoint("Receiver::Send Galois Keys");

    co_await send_encoded_chunks(FX, chl);

    vector<Ciphertext> decoded_in_he;
    co_await recv_decoded_chunks(decoded_in_he, chl);

    co_await chl.recvResize(last_layer_per_bin);
    auto len = 0;
    for (size_t i = 0; i < mNumSlots; i++)
    {
      len += last_layer_per_bin[i];
    }
    oc::BitVector occupy_indicator_flat(len);
    co_await chl.recv(occupy_indicator_flat);
    auto offset = 0;
    mLayerToBins.clear();
    mLayerToBins.resize(decoded_in_he.size());
    for (size_t i = 0; i < mNumSlots; i++)
    {
      for (uint32_t layer = 0; layer < last_layer_per_bin[i]; ++layer)
      {
        if (occupy_indicator_flat[offset + layer])
        {
          mLayerToBins[layer].push_back(static_cast<uint32_t>(i));
        }
      }
      offset += last_layer_per_bin[i];
    }

    vector<uint64_t> dec_results;
    decrypt(decoded_in_he, dec_results);

    if (!mSharedOutput)
    {
      results.resize(mNsender);
      for (size_t i = 0; i < mNsender; i++)
      {
        results[i] = (dec_results[i] == mIndicatorStr) ? 1 : 0;
      }
    }
    else
    {
      // GMW with dec_results
      u64 keyBitLength = 40 + oc::log2ceil(mNsender);
      u64 keyByteLength = oc::divCeil(keyBitLength, 8);

      oc::Matrix<u8> gmwin;

      gmwin.resize(mNsender, keyByteLength, oc::AllocType::Uninitialized);
      for (size_t i = 0; i < mNsender; i++)
      {
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

  Proto RpmtReceiver::send_encoded_chunks(const std::vector<oc::block> &X,
                                          Socket &chl)
  {
    mIndicatorStr =
        seal::util::barrett_reduce_64(mPrng.get<uint64_t>(), mModulus);
    vector<uint64_t> val(mN, mIndicatorStr);

    PrimeFieldOkvs okvs;
    // okvs.setTimer(getTimer());
    okvs.init(X.size(), mM, mW, mModulus);

    vector<uint64_t> encoded(mM);

    if (!okvs.encode(X, val, encoded))
      throw RTE_LOC;

    vector<uint64_t> encoded_spacing(mM);
    for (uint32_t i = 0; i < encoded.size(); i++)
    {
      uint32_t q = i / mHalfSlots;
      uint32_t r = i % mHalfSlots;
      encoded_spacing[i] = encoded[r * mNumHalfBatch + q];
    }
    encoded = encoded_spacing;

    setTimePoint("Receiver::OKVS Encoding");

    vector<uint64_t> plainVec(mNumSlots);
    Plaintext ptxt;
    stringstream ctxtstream;
    const uint64_t sentCtxts = totalEncodedCipherCount(mNumBatch, mW, mWrap);

    const uint32_t chunkSize = encodedBatchChunk();
    for (uint32_t jBegin = 0; jBegin < mNumBatch; jBegin += chunkSize)
    {
      const uint32_t jEnd = std::min<uint32_t>(mNumBatch, jBegin + chunkSize);
      ctxtstream.clear();
      ctxtstream.str("");

      for (uint32_t j = jBegin; j < jEnd; ++j)
      {
        std::copy_n(&encoded[j * mNumSlots], mNumSlots, plainVec.data());
        mBatchEncoder->encode(plainVec, ptxt);
        mEncryptor->encrypt_symmetric(ptxt).save(ctxtstream);
      }

      auto payload = ctxtstream.str();
      co_await chl.send(move(payload));
    }

    cout << "Receiver sends " << sentCtxts
         << " Ctxts of OKVS Encoding" << endl;

    setTimePoint("Receiver::Encryption & Send");
  }

  Proto RpmtReceiver::recv_decoded_chunks(
      std::vector<seal::Ciphertext> &decoded_in_he,
      Socket &chl)
  {
    uint32_t decoded_he_size;
    co_await chl.recv(decoded_he_size);

    decoded_in_he.resize(decoded_he_size);
    SEALContext context = *mContext;
    string recvstring;
    stringstream recvstream;

    const uint32_t chunkSize = decodedLayerChunk();
    for (uint32_t layerBegin = 0; layerBegin < decoded_he_size; layerBegin += chunkSize)
    {
      const uint32_t layerEnd =
          std::min<uint32_t>(decoded_he_size, layerBegin + chunkSize);
      co_await chl.recvResize(recvstring);

      recvstream.clear();
      recvstream.str(recvstring);

      for (uint32_t i = layerBegin; i < layerEnd; ++i)
      {
        decoded_in_he[i].unsafe_load(context, recvstream);
      }
    }

    setTimePoint("Receiver::Recv back and Serialize");
  }

  void RpmtReceiver::decrypt(const std::vector<seal::Ciphertext> &decoded_in_he,
                              std::vector<uint64_t> &dec_results)
  {
    const size_t L = decoded_in_he.size();
    vector<Plaintext> ptxts(L);

    cout << "Noise Budget: "
         << mDecryptor->invariant_noise_budget(decoded_in_he[0]) << endl;

    for (size_t i = 0; i < L; i++)
    {
      mDecryptor->decrypt(decoded_in_he[i], ptxts[i]);
    }

    dec_results.resize(mNsender);
    vector<uint64_t> decodeVec(mNumSlots);

    auto idx = 0;
    for (size_t i = 0; i < L; i++)
    {
      mBatchEncoder->decode(ptxts[i], decodeVec);

      for (auto bin : mLayerToBins[i])
      {
        dec_results[idx++] = decodeVec[bin];
      }
    }
    setTimePoint("Receiver::Decrypt");
  }
} // namespace rlweOkvs
