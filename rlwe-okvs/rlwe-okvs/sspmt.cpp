#include "sspmt.h"

#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <utility>

#include "../GMW/Gmw.h"
#include "okvs.h"
#include "seal/util/numth.h"
#include "seal/util/uintarithsmallmod.h"

#include "macoro/when_all.h"

using namespace std;
using namespace seal;
using namespace oc;
using namespace volePSI;

namespace rlweOkvs {
namespace {
uint32_t chunkEnvOrDefault(const char *name, uint32_t fallback) {
  if (const char *value = std::getenv(name)) {
    auto parsed = static_cast<uint32_t>(std::strtoul(value, nullptr, 10));
    if (parsed > 0) {
      return parsed;
    }
  }
  return fallback;
}

uint32_t encodedBatchChunk() {
  return chunkEnvOrDefault("RLWE_OKVS_ENCODED_BATCH_CHUNK", 8);
}

uint32_t decodedLayerChunk() {
  return chunkEnvOrDefault("RLWE_OKVS_DECODED_LAYER_CHUNK", 8);
}

uint32_t activeBatchCountForWrap(uint32_t k, uint32_t numBatch,
                                 uint32_t width) {
  const int64_t count =
      static_cast<int64_t>(width) + (1 - static_cast<int64_t>(k)) * numBatch - 1;
  if (count <= 0) {
    return 0;
  }
  return std::min<uint32_t>(numBatch, static_cast<uint32_t>(count));
}

uint32_t activeWrapCountForBatch(uint32_t j, uint32_t numBatch, uint32_t width,
                                 uint32_t maxWrap) {
  uint32_t wraps = 0;
  while (wraps < maxWrap && j < activeBatchCountForWrap(wraps, numBatch, width)) {
    ++wraps;
  }
  return wraps;
}

// log P[Bin(n, p) >= k], summed term by term in log space. Exact up to the
// truncation of terms more than e^-80 below the largest one (a relative error
// below 2e-35), which is far tighter than the 2^-70-ish levels we test.
double logBinomUpperTail(uint64_t n, double logP, double logQ, double lgn1,
                         uint64_t k, double mean) {
  if (k == 0) return 0.0;
  if (k > n) return -std::numeric_limits<double>::infinity();

  double maxLog = -std::numeric_limits<double>::infinity();
  double sum = 0.0;
  for (uint64_t j = k; j <= n; ++j) {
    const double l = lgn1 - std::lgamma((double)j + 1.0) -
                     std::lgamma((double)(n - j) + 1.0) + (double)j * logP +
                     (double)(n - j) * logQ;
    if (l > maxLog) {
      sum = (maxLog == -std::numeric_limits<double>::infinity())
                ? 1.0
                : sum * std::exp(maxLog - l) + 1.0;
      maxLog = l;
    } else {
      sum += std::exp(l - maxLog);
    }
    if ((double)j > mean && l < maxLog - 80.0) break;
  }
  return maxLog + std::log(sum);
}

uint64_t totalEncodedCipherCount(uint32_t numBatch, uint32_t width,
                                 uint32_t maxWrap) {
  uint64_t total = 0;
  for (uint32_t k = 0; k < maxWrap; ++k) {
    total += activeBatchCountForWrap(k, numBatch, width);
  }
  return total;
}
}

uint32_t SspmtSender::sequenceLayers(const std::vector<uint32_t> &itemBin,
                                     const std::vector<uint32_t> &itemBlock,
                                     uint32_t numSlots, uint32_t spanBlocks,
                                     std::vector<uint32_t> &itemToLayer,
                                     std::vector<uint32_t> &layerMinBlock,
                                     std::vector<uint32_t> &layerMaxBlock) {
  const uint32_t n = static_cast<uint32_t>(itemBin.size());
  itemToLayer.assign(n, UINT32_MAX);
  layerMinBlock.clear();
  layerMaxBlock.clear();
  if (n == 0) {
    return 0;
  }

  uint32_t max_block = 0;
  for (uint32_t i = 0; i < n; ++i) {
    max_block = std::max(max_block, itemBlock[i]);
  }

  std::vector<std::vector<uint32_t>> block_items(max_block + 1);
  for (uint32_t i = 0; i < n; ++i) {
    block_items[itemBlock[i]].push_back(i);
  }

  // Blocks are processed left to right, so a layer's min block is fixed at
  // creation and span admissibility is purely "min_block + spanBlocks > blk".
  // Layers are created with nondecreasing min block; the admissible ones form
  // a suffix [firstActive, end), and scanning it in order is exactly
  // "smallest anchor first" (earliest-expiring first).
  std::vector<std::vector<uint8_t>> used_bins;
  uint32_t firstActive = 0;

  for (uint32_t blk = 0; blk <= max_block; ++blk) {
    while (firstActive < layerMinBlock.size() &&
           layerMinBlock[firstActive] + spanBlocks <= blk) {
      ++firstActive;
    }
    for (uint32_t idx : block_items[blk]) {
      const uint32_t r = itemBin[idx];
      bool placed = false;
      for (uint32_t li = firstActive; li < used_bins.size(); ++li) {
        if (used_bins[li][r]) {
          continue;
        }
        used_bins[li][r] = 1;
        layerMaxBlock[li] = blk;
        itemToLayer[idx] = li;
        placed = true;
        break;
      }
      if (!placed) {
        used_bins.emplace_back(numSlots, 0);
        used_bins.back()[r] = 1;
        layerMinBlock.push_back(blk);
        layerMaxBlock.push_back(blk);
        itemToLayer[idx] = static_cast<uint32_t>(used_bins.size()) - 1;
      }
    }
  }
  return static_cast<uint32_t>(used_bins.size());
}

uint32_t SspmtSender::sequenceLayersOptimal(
    const std::vector<uint32_t> &itemBin,
    const std::vector<uint32_t> &itemBlock, uint32_t numSlots,
    uint32_t numBlocks, uint32_t spanBlocks,
    std::vector<uint32_t> &itemToLayer, std::vector<uint32_t> &layerMinBlock,
    std::vector<uint32_t> &layerMaxBlock) {
  const uint32_t n = static_cast<uint32_t>(itemBin.size());
  itemToLayer.assign(n, UINT32_MAX);
  layerMinBlock.clear();
  layerMaxBlock.clear();
  if (n == 0) {
    return 0;
  }
  const uint32_t b = numBlocks;
  const uint32_t z = std::max<uint32_t>(spanBlocks, 1);

  std::vector<std::vector<uint32_t>> block_items(b);
  for (uint32_t i = 0; i < n; ++i) {
    block_items[itemBlock[i]].push_back(i);
  }

  // --- Phase 1: how many windows, and where ---------------------------------
  // cnt[x][bin] counts the items of that bin with block in [x, y] as y sweeps,
  // and Mrun[x] = M([x, y]). pref[s] is the number of windows placed at starts
  // <= s (final for s < y, since windows are only ever added at the current
  // y). The supply reaching an interval [x, y] is total - pref[x-z], so the
  // shortfall over all x is max_x (Mrun[x] + pref[x-z]) - total.
  std::vector<std::vector<uint32_t>> cnt(b, std::vector<uint32_t>(numSlots, 0));
  std::vector<uint32_t> Mrun(b, 0);
  std::vector<uint64_t> pref(b, 0);
  std::vector<uint32_t> anchors;
  uint64_t total = 0;

  for (uint32_t y = 0; y < b; ++y) {
    for (uint32_t idx : block_items[y]) {
      const uint32_t bin = itemBin[idx];
      for (uint32_t x = 0; x <= y; ++x) {
        const uint32_t c = ++cnt[x][bin];
        if (c > Mrun[x]) {
          Mrun[x] = c;
        }
      }
    }

    uint64_t need = 0;
    for (uint32_t x = 0; x <= y; ++x) {
      const uint64_t supplyBelow = (x >= z) ? pref[x - z] : 0;
      need = std::max(need, (uint64_t)Mrun[x] + supplyBelow);
    }
    if (need > total) {
      anchors.insert(anchors.end(), need - total, y);
      total = need;
    }
    pref[y] = total;
  }

  // --- Phase 2: assign, leftmost admissible window first ---------------------
  // Anchors are nondecreasing and each bin's items arrive in block order, so a
  // single forward pointer per bin suffices: windows before it are either used
  // by this bin or can no longer cover its remaining items.
  const uint32_t L = static_cast<uint32_t>(anchors.size());
  std::vector<uint32_t> ptr(numSlots, 0);
  std::vector<uint32_t> lmin(L, UINT32_MAX), lmax(L, 0);

  for (uint32_t y = 0; y < b; ++y) {
    for (uint32_t idx : block_items[y]) {
      const uint32_t bin = itemBin[idx];
      uint32_t &p = ptr[bin];
      const uint32_t lo = (y >= z - 1) ? (y - z + 1) : 0;
      while (p < L && anchors[p] < lo) {
        ++p;
      }
      if (p == L || anchors[p] > y) {
        // Only reachable if a Hall constraint was violated, which phase 1
        // rules out.
        throw RTE_LOC;
      }
      itemToLayer[idx] = p;
      lmin[p] = std::min(lmin[p], y);
      lmax[p] = std::max(lmax[p], y);
      ++p;
    }
  }

  // --- Phase 3: drop the windows no item landed in --------------------------
  std::vector<uint32_t> remap(L, UINT32_MAX);
  uint32_t kept = 0;
  for (uint32_t l = 0; l < L; ++l) {
    if (lmin[l] != UINT32_MAX) {
      remap[l] = kept++;
      layerMinBlock.push_back(lmin[l]);
      layerMaxBlock.push_back(lmax[l]);
    }
  }
  for (uint32_t i = 0; i < n; ++i) {
    itemToLayer[i] = remap[itemToLayer[i]];
  }
  return kept;
}

u32 sspmtParams::resolveLayerBudget(u64 n) const {
  if (layerBudget) {
    return layerBudget;
  }
  if (n == 0) {
    return 0;
  }

  const uint32_t numSlots = heNumSlots;
  const uint64_t m = roundUpTo(bandExpansion * n, numSlots);
  const uint32_t b = (uint32_t)(m / numSlots);
  const uint64_t positionRange = m - bandWidth + 1;
  const uint32_t lambda = layerBudgetLambda;
  const uint32_t W = std::max<uint32_t>(span_blocks, 1);
  if (b == 0) {
    return 0;
  }

  // One share of the 2^-lambda budget per (bin, start position, length)
  // Hall constraint.
  const double logEps = -(double)lambda * std::log(2.0) -
                        std::log((double)numSlots) -
                        2.0 * std::log((double)b);
  const double lgn1 = std::lgamma((double)n + 1.0);

  uint64_t k = 0;  // D(g); nondecreasing in g, so the walk is amortized
  uint64_t budget = 0;
  for (uint32_t g = 1; g <= b; ++g) {
    const double p = (double)g / (double)positionRange;
    const double logP = std::log((double)g) - std::log((double)positionRange);
    const double logQ = std::log1p(-p);
    const double mean = (double)n * p;

    while (logBinomUpperTail(n, logP, logQ, lgn1, k, mean) > logEps) {
      ++k;
    }

    // A layer spans W blocks, so it can serve an item at block beta iff its
    // start lies in [beta-W+1, beta]: the starts that reach a length-g
    // interval number g+W-1, and the starts covering all b blocks number
    // b+W-1. Hence rho >= (D(g)+1)/(g+W-1) for every g, and the window count
    // is ceil(rho * (b+W-1)); kept in integers to avoid rounding surprises.
    const uint64_t val =
        divCeil((k + 1) * (uint64_t)(b + W - 1), (uint64_t)(g + W - 1));
    budget = std::max(budget, val);
  }
  return (uint32_t)budget;
}

void SspmtSender::sequencing(const std::vector<uint32_t> &start_pos_spacing) {
  std::vector<uint32_t> item_binidx(mN);
  for (uint32_t i = 0; i < mN; ++i) {
    uint32_t pos = start_pos_spacing[i];
    item_binidx[i] = pos % mNumSlots;
    mItemToBlockIdx[i] = pos / mNumSlots;
  }

  mNumRealLayers = sequenceLayers(item_binidx, mItemToBlockIdx, mNumSlots,
                                  mSpanBlocks, mItemToLayerIdx, mLayerMinBlock,
                                  mLayerMaxBlock);

  // The realized layer count is a function of Y (its collision structure), so
  // it is padded up to the public budget before anything leaves this party.
  //
  // The budget bounds the *minimum* layer count, so overshooting it with the
  // greedy is not yet the abort event: fall back to the optimal sequencing
  // first. (The greedy has matched the optimum on every instance we have
  // measured, so this path is not expected to run -- but the 2^-lambda abort
  // guarantee should not rest on that.)
  if (mNumRealLayers > mLayerBudget) {
    mNumRealLayers = sequenceLayersOptimal(
        item_binidx, mItemToBlockIdx, mNumSlots, mNumBatch, mSpanBlocks,
        mItemToLayerIdx, mLayerMinBlock, mLayerMaxBlock);
  }
  if (mNumRealLayers > mLayerBudget) {
    std::cout << "sequencing: " << mNumRealLayers << " layers exceed the "
              << mLayerBudget << " layer budget" << std::endl;
    throw RTE_LOC;
  }

  // The layout. Padding layers hold no item; their anchor is bookkeeping only.
  // They are left where the sequencer put them, in a suffix: the sequencer
  // front-loads each bin, so early layers are denser, but that is a prior the
  // receiver can compute from the public parameters alone -- its view is
  // uniform shares, with nothing in it that correlates with occupancy, so
  // there is nothing for the prior to be applied to. (Where slot coordinates
  // do become visible, as in the PSU transfer, permuting the layers would not
  // be enough either; that needs the permute+share step noted in pso.h.)
  mSlotToItem.assign(static_cast<size_t>(mNumLayers) * mNumSlots, UINT32_MAX);
  for (uint32_t i = 0; i < mN; ++i) {
    mSlotToItem[static_cast<size_t>(mItemToLayerIdx[i]) * mNumSlots +
                item_binidx[i]] = i;
  }

  mLayerMinBlock.resize(mNumLayers);
  mLayerMaxBlock.resize(mNumLayers);
  mLayerIsPadding.assign(mNumLayers, 0);
  for (uint32_t l = mNumRealLayers; l < mNumLayers; ++l) {
    mLayerMinBlock[l] = 0;
    mLayerMaxBlock[l] = 0;
    mLayerIsPadding[l] = 1;
  }

  // Every slot of the L x H layout -- occupied, empty or padding -- gets a
  // mask and takes part in the equality, so the occupancy is never revealed.
  // An empty slot of a real layer decodes to 0, so the sender's mask r and the
  // receiver's r - indicator differ (the indicator is nonzero) and the slot
  // shares a 0 deterministically; a padding slot decodes to a random value and
  // only collides with the indicator with probability 2^-log(p) per slot.
  mMasks.resize(static_cast<size_t>(mNumLayers) * mNumSlots);
  ptxts_mask.resize(mNumLayers);
  vector<uint64_t> raw_masks(mNumSlots);
  for (size_t lay = 0; lay < mNumLayers; lay++) {
    mPrng.get<uint64_t>(raw_masks);
    for (uint32_t bin = 0; bin < mNumSlots; bin++) {
      raw_masks[bin] = seal::util::barrett_reduce_64(raw_masks[bin], mModulus);
      mMasks[lay * mNumSlots + bin] = raw_masks[bin];
    }
    mBatchEncoder->encode(raw_masks, ptxts_mask[lay]);
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

  // The layout is fixed by the public parameters, not by the sequencing --
  // the sender pads up to it. Setting the layer count here rather than in
  // sequencing() is what lets setup() size the GMW before any input is known.
  mLayerBudget = ssParams.resolveLayerBudget(mN);
  mNumLayers = mLayerBudget;

  mItemToBlockIdx.resize(mN);
  mItemToLayerIdx.resize(mN);

  parms.set_coeff_modulus(
      CoeffModulus::Create(mNumSlots, ssParams.heCoeffModulus));
  mModulus =
      seal::util::get_primes(2 * mNumSlots, ssParams.hePlainModulusBits, 4)[3];
  parms.set_plain_modulus(mModulus);

  mContext = make_shared<SEALContext>(parms, true, sec_level_type::none);
  mBatchEncoder = make_unique<BatchEncoder>(*mContext);
  mEvaluator = make_unique<Evaluator>(*mContext);
};

void SspmtSender::initGmw() {
  const u64 keyBitLength = 40 + oc::log2ceil(mN);
  auto cir = isZeroCircuit(keyBitLength);
  mGmw.setTimer(getTimer());
  mGmw.init(getLayoutSize(), cir, 1, mOTeBatchSize, 0, mPrng.get());
}

Proto SspmtSender::setup(Socket &chl) {
  if (mSetupDone) {
    co_return;
  }
  initGmw();
  co_await mGmw.generateTriple(chl);
  mSetupDone = true;
  setTimePoint("Sender::Setup (GMW triples)");
}

Proto SspmtSender::run(const std::vector<oc::block> &Y, oc::BitVector &results,
                       Socket &chl) {
  // The triples are input independent; if the offline phase was skipped they
  // are generated here instead, before any input-dependent work starts.
  co_await setup(chl);

  preprocess(Y);
  vector<vector<Ciphertext>> encoded_in_he(mNumBatch);
  co_await recv_encoded_chunks(encoded_in_he, chl);
  co_await send_decoded_chunks(encoded_in_he, chl);

  const u64 nInst = getLayoutSize();
  const u64 keyBitLength = 40 + oc::log2ceil(mN);
  const u64 keyByteLength = oc::divCeil(keyBitLength, 8);

  oc::Matrix<u8> gmwin(nInst, keyByteLength, oc::AllocType::Uninitialized);
  for (size_t i = 0; i < nInst; i++) {
    memcpy(&gmwin(i, 0), &mMasks[i], keyByteLength);
  }
  mGmw.setInput(0, gmwin);
  co_await mGmw.run(chl);  // triples already generated -> online only

  auto rr = mGmw.getOutputView(0);
  results.resize(nInst);
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

  struct BinMeta {
    uint32_t bin;
    const uint64_t *band_ptr;
    uint32_t start_blk;
  };

  using Contrib = std::pair<uint32_t, uint64_t>;
  constexpr uint32_t B_CHUNK = 128;

  const uint64_t *bands = bands_flat.data();

  ptxts_diags.resize(mNumLayers);
  for (size_t i = 0; i < mNumLayers; ++i) {
    ptxts_diags[i].resize(static_cast<size_t>(mNumBatch) * mWrap);
  }

  std::vector<uint64_t> plainVec(mNumSlots, 0);
  std::vector<BinMeta> layer_meta;
  layer_meta.reserve(mNumSlots);
  std::vector<uint32_t> row_counts(B_CHUNK);
  std::vector<uint32_t> row_offsets(B_CHUNK + 1);
  std::vector<uint32_t> write_ptr(B_CHUNK + 1);
  std::vector<Contrib> flat_contribs;

  for (uint32_t i = 0; i < mNumLayers; ++i) {
    // Padding layers get their (single) random multiplier in
    // encrypted_decode(); nothing to prepare here.
    if (mLayerIsPadding[i]) {
      continue;
    }

    uint32_t Bmin = mLayerMinBlock[i];
    uint32_t Bmax = mLayerMaxBlock[i] + (mW - 1);

    layer_meta.clear();

    for (uint32_t bin = 0; bin < mNumSlots; ++bin) {
      uint32_t item = mSlotToItem[static_cast<size_t>(i) * mNumSlots + bin];
      if (item == UINT32_MAX) {
        continue;
      }

      layer_meta.push_back(
          BinMeta{bin, bands + static_cast<size_t>(item) * mW,
                  mItemToBlockIdx[item]});
    }

    for (uint32_t chunk_begin = Bmin; chunk_begin <= Bmax;
         chunk_begin += B_CHUNK) {
      const uint32_t chunk_end =
          std::min<uint32_t>(Bmax, chunk_begin + B_CHUNK - 1);
      const uint32_t chunk_len = chunk_end - chunk_begin + 1;

      std::fill_n(row_counts.begin(), chunk_len, 0);

      for (const auto &bm : layer_meta) {
        const uint32_t start = bm.start_blk;
        const uint32_t end = bm.start_blk + mW - 1;
        const uint32_t lo = std::max(start, chunk_begin);
        const uint32_t hi = std::min(end, chunk_end);

        if (lo > hi) {
          continue;
        }

        for (uint32_t B = lo; B <= hi; ++B) {
          ++row_counts[B - chunk_begin];
        }
      }

      row_offsets[0] = 0;
      for (uint32_t row = 0; row < chunk_len; ++row) {
        row_offsets[row + 1] = row_offsets[row] + row_counts[row];
      }

      flat_contribs.resize(row_offsets[chunk_len]);
      std::copy_n(row_offsets.begin(), chunk_len + 1, write_ptr.begin());

      for (const auto &bm : layer_meta) {
        const uint32_t start = bm.start_blk;
        const uint32_t end = bm.start_blk + mW - 1;
        const uint32_t lo = std::max(start, chunk_begin);
        const uint32_t hi = std::min(end, chunk_end);

        if (lo > hi) {
          continue;
        }

        for (uint32_t B = lo; B <= hi; ++B) {
          const uint32_t row = B - chunk_begin;
          const uint32_t w = B - start;
          flat_contribs[write_ptr[row]++] = Contrib{bm.bin, bm.band_ptr[w]};
        }
      }

      for (uint32_t row = 0; row < chunk_len; ++row) {
        const uint32_t B = chunk_begin + row;
        const uint32_t k = B / mNumBatch;
        const uint32_t j = B % mNumBatch;
        if (k >= mWrap) {
          continue;
        }

        const size_t outIdx = static_cast<size_t>(j) * mWrap + k;
        const uint32_t begin = row_offsets[row];
        const uint32_t end = row_offsets[row + 1];

        for (uint32_t t = begin; t < end; ++t) {
          const auto &cv = flat_contribs[t];
          plainVec[cv.first] = cv.second;
        }

        mBatchEncoder->encode(plainVec, ptxts_diags[i][outIdx]);

        for (uint32_t t = begin; t < end; ++t) {
          const auto &cv = flat_contribs[t];
          plainVec[cv.first] = 0;
        }
      }
    }
  }
  setTimePoint("Sender::HE Encode");
}

void SspmtSender::encrypted_decode(
    const std::vector<std::vector<seal::Ciphertext>> &encoded_in_he,
    std::vector<seal::Ciphertext> &decoded_in_he, uint32_t layerBegin,
    uint32_t layerEnd) {
  decoded_in_he.resize(layerEnd - layerBegin);

  std::vector<uint64_t> padVec(mNumSlots);
  Plaintext padPtxt;

  for (uint32_t i = layerBegin; i < layerEnd; ++i) {
    bool initialized = false;
    Ciphertext tmp;
    Ciphertext &out = decoded_in_he[i - layerBegin];

    // A padding layer's ciphertext only has to decrypt to something the
    // receiver cannot predict: one received ciphertext times a dense random
    // plaintext, plus the mask. A single multiply instead of a real layer's
    // w + span chain -- the padding cost is compute-negligible. Its noise is
    // therefore smaller than a real layer's; noise-level indistinguishability
    // (which pre-flooding is broken anyway, since a real layer's noise scales
    // with its occupancy) is the job of the upcoming noise-flooding step.
    if (mLayerIsPadding[i]) {
      mPrng.get<uint64_t>(padVec);
      for (uint32_t bin = 0; bin < mNumSlots; ++bin) {
        padVec[bin] = seal::util::barrett_reduce_64(padVec[bin], mModulus);
      }
      mBatchEncoder->encode(padVec, padPtxt);
      mEvaluator->multiply_plain(encoded_in_he[0][0], padPtxt, out);
      mEvaluator->add_plain_inplace(out, ptxts_mask[i]);
      mEvaluator->mod_switch_to_next_inplace(out);
      continue;
    }

    const uint32_t Bmin = mLayerMinBlock[i];
    const uint32_t Bmax = mLayerMaxBlock[i] + (mW - 1);

    for (uint32_t j = 0; j < mNumBatch; ++j) {
      if (Bmax < j) {
        continue;
      }

      uint32_t k_begin = 0;
      if (Bmin > j) {
        k_begin = (Bmin - j + mNumBatch - 1) / mNumBatch;
      }

      if (k_begin >= mWrap) {
        continue;
      }

      uint32_t k_end = (Bmax - j) / mNumBatch;
      if (k_end >= mWrap) {
        k_end = static_cast<uint32_t>(mWrap - 1);
      }

      if (k_begin > k_end) {
        continue;
      }

      const uint32_t activeWraps =
          activeWrapCountForBatch(j, mNumBatch, mW, mWrap);
      if (k_begin >= activeWraps) {
        continue;
      }
      k_end = std::min<uint32_t>(k_end, activeWraps - 1);

      uint32_t k = k_begin;

      if (!initialized) {
        const size_t diag_idx = static_cast<size_t>(j) * mWrap + k;
        mEvaluator->multiply_plain(encoded_in_he[j][k],
                                   ptxts_diags[i][diag_idx], out);
        initialized = true;
        ++k;
      }

      for (; k <= k_end; ++k) {
        const size_t diag_idx = static_cast<size_t>(j) * mWrap + k;
        mEvaluator->multiply_plain(encoded_in_he[j][k],
                                   ptxts_diags[i][diag_idx], tmp);
        mEvaluator->add_inplace(out, tmp);
      }
    }

    mEvaluator->add_plain_inplace(out, ptxts_mask[i]);
    mEvaluator->mod_switch_to_next_inplace(out);
  }
}

Proto SspmtSender::recv_encoded_chunks(
    std::vector<std::vector<seal::Ciphertext>> &encoded_in_he, Socket &chl) {
  SEALContext context = *mContext;
  size_t recvBytes = 0;
  uint64_t recvCtxts = 0;
  string recvstring;
  stringstream recvstream;

  const uint32_t chunkSize = encodedBatchChunk();
  for (uint32_t jBegin = 0; jBegin < mNumBatch; jBegin += chunkSize) {
    const uint32_t jEnd = std::min<uint32_t>(mNumBatch, jBegin + chunkSize);
    co_await chl.recvResize(recvstring);
    recvBytes += recvstring.size();

    recvstream.clear();
    recvstream.str(recvstring);

    for (uint32_t j = jBegin; j < jEnd; ++j) {
      encoded_in_he[j].resize(mWrap);
      const uint32_t activeWraps =
          activeWrapCountForBatch(j, mNumBatch, mW, mWrap);
      for (uint32_t k = 0; k < activeWraps; ++k) {
        encoded_in_he[j][k].unsafe_load(context, recvstream);
        ++recvCtxts;
      }
    }
  }

  cout << "Sender receives " << recvCtxts << " Ctxts of OKVS Encoding, "
       << recvBytes << " Bytes" << endl;

  setTimePoint("Sender::Recv ctxts & Serialize");
}

Proto SspmtSender::send_decoded_chunks(
    const std::vector<std::vector<seal::Ciphertext>> &encoded_in_he,
    Socket &chl) {
  // The layer count is the public budget, known to both parties from the
  // parameters, so it is never transmitted: the realized count depends on Y.
  size_t sentBytes = 0;
  std::vector<Ciphertext> decoded_chunk;
  stringstream sendstream;
  const uint32_t chunkSize = decodedLayerChunk();

  for (uint32_t layerBegin = 0; layerBegin < mNumLayers;
       layerBegin += chunkSize) {
    const uint32_t layerEnd =
        std::min<uint32_t>(mNumLayers, layerBegin + chunkSize);

    encrypted_decode(encoded_in_he, decoded_chunk, layerBegin, layerEnd);

    sendstream.clear();
    sendstream.str("");
    for (size_t i = 0; i < decoded_chunk.size(); ++i) {
      decoded_chunk[i].save(sendstream);
    }

    auto payload = sendstream.str();
    sentBytes += payload.size();
    co_await chl.send(move(payload));
  }

  cout << "Sender sends " << mNumLayers << " decoded ctxts, " << sentBytes
       << " Bytes" << endl;

  setTimePoint("Sender::Encrypted OKVS Decoding & Send Back");
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

  // Same public formula as the sender's, over the sender's set size: both
  // parties must agree on the layout size without communicating it.
  mLayerBudget = ssParams.resolveLayerBudget(mNsender);

  parms.set_coeff_modulus(
      CoeffModulus::Create(mNumSlots, ssParams.heCoeffModulus));
  mModulus =
      seal::util::get_primes(2 * mNumSlots, ssParams.hePlainModulusBits, 4)[3];
  parms.set_plain_modulus(mModulus);

  mContext = make_shared<SEALContext>(parms, true, sec_level_type::none);
  KeyGenerator keygen(*mContext);
  SecretKey secret_key = keygen.secret_key();
  PublicKey public_key;
  keygen.create_public_key(public_key);

  mEncryptor = make_unique<Encryptor>(*mContext, public_key);
  mEncryptor->set_secret_key(secret_key);

  mBatchEncoder = make_unique<BatchEncoder>(*mContext);
  mDecryptor = make_unique<Decryptor>(*mContext, secret_key);
};

void SspmtReceiver::initGmw() {
  const u64 keyBitLength = 40 + oc::log2ceil(mNsender);
  auto cir = isZeroCircuit(keyBitLength);
  mGmw.setTimer(getTimer());
  mGmw.init(getLayoutSize(), cir, 1, mOTeBatchSize, 1, mPrng.get());
}

Proto SspmtReceiver::setup(Socket &chl) {
  if (mSetupDone) {
    co_return;
  }
  initGmw();
  co_await mGmw.generateTriple(chl);
  mSetupDone = true;
  setTimePoint("Receiver::Setup (GMW triples)");
}

Proto SspmtReceiver::run(const std::vector<oc::block> &X,
                         oc::BitVector &results, Socket &chl) {
  co_await setup(chl);

  co_await send_encoded_chunks(X, chl);

  // The equality runs over the whole L x H rectangle -- no occupancy is
  // received -- so slot (layer i, bin b) matches iff the Y-item sitting there
  // (if any) is in X; empty slots decode to 0 and never match.
  vector<Ciphertext> decoded_in_he;
  co_await recv_decoded_chunks(decoded_in_he, chl);

  const size_t nInst = getLayoutSize();
  const u64 keyBitLength = 40 + oc::log2ceil(mNsender);
  const u64 keyByteLength = oc::divCeil(keyBitLength, 8);

  // How much room is left for the noise flooding that re-randomization will
  // need: the returned ciphertexts must still decrypt after e_flood is added,
  // so the smallest budget here is the ceiling on log2(e_flood / e_compute).
  {
    int minBudget = INT32_MAX, maxBudget = 0;
    for (size_t i = 0; i < mLayerBudget; ++i) {
      const int budget = mDecryptor->invariant_noise_budget(decoded_in_he[i]);
      minBudget = std::min(minBudget, budget);
      maxBudget = std::max(maxBudget, budget);
    }
    cout << "Return noise budget: " << minBudget << " - " << maxBudget
         << " bits" << endl;
  }

  oc::Matrix<u8> gmwin(nInst, keyByteLength, oc::AllocType::Uninitialized);
  vector<uint64_t> decodeVec(mNumSlots);
  Plaintext ptxt;
  for (size_t i = 0; i < mLayerBudget; ++i) {
    mDecryptor->decrypt(decoded_in_he[i], ptxt);
    mBatchEncoder->decode(ptxt, decodeVec);
    for (uint32_t bin = 0; bin < mNumSlots; ++bin) {
      auto tmp = util::sub_uint_mod(decodeVec[bin], mIndicatorStr, mModulus);
      memcpy(&gmwin(i * mNumSlots + bin, 0), &tmp, keyByteLength);
    }
  }
  setTimePoint("Receiver::Decrypt");

  mGmw.setInput(0, gmwin);
  co_await mGmw.run(chl);  // triples already generated -> online only

  auto rr = mGmw.getOutputView(0);
  results.resize(nInst);
  std::copy(rr.begin(), rr.end(), results.data());
  setTimePoint("Receiver::Online GMW");
}

Proto SspmtReceiver::send_encoded_chunks(const std::vector<oc::block> &X,
                                         Socket &chl) {
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
  stringstream ctxtstream;
  const uint64_t sentCtxts = totalEncodedCipherCount(mNumBatch, mW, mWrap);

  const uint32_t chunkSize = encodedBatchChunk();
  for (uint32_t jBegin = 0; jBegin < mNumBatch; jBegin += chunkSize) {
    const uint32_t jEnd = std::min<uint32_t>(mNumBatch, jBegin + chunkSize);
    ctxtstream.clear();
    ctxtstream.str("");

    for (uint32_t j = jBegin; j < jEnd; ++j) {
      const uint32_t activeWraps =
          activeWrapCountForBatch(j, mNumBatch, mW, mWrap);
      for (uint32_t k = 0; k < activeWraps; ++k) {
        if (j == mNumBatch - 1 && k > 0) {
          std::copy_n(&encoded[j * mNumSlots + k], mNumSlots - k,
                      plainVec.data());
          std::copy_n(&encoded[0], k, plainVec.data() + mNumSlots - k);
        } else {
          std::copy_n(&encoded[j * mNumSlots + k], mNumSlots, plainVec.data());
        }
        mBatchEncoder->encode(plainVec, ptxt);
        mEncryptor->encrypt_symmetric(ptxt).save(ctxtstream);
      }
    }

    auto payload = ctxtstream.str();
    co_await chl.send(move(payload));
  }

  cout << "Receiver sends " << sentCtxts << " Ctxts of OKVS Encoding" << endl;

  setTimePoint("Receiver::Encryption & Send");
}

Proto SspmtReceiver::recv_decoded_chunks(
    std::vector<seal::Ciphertext> &decoded_in_he, Socket &chl) {
  const uint32_t numLayers = mLayerBudget;
  decoded_in_he.resize(numLayers);
  SEALContext context = *mContext;
  string recvstring;
  stringstream recvstream;

  const uint32_t chunkSize = decodedLayerChunk();
  for (uint32_t layerBegin = 0; layerBegin < numLayers;
       layerBegin += chunkSize) {
    const uint32_t layerEnd =
        std::min<uint32_t>(numLayers, layerBegin + chunkSize);
    co_await chl.recvResize(recvstring);

    recvstream.clear();
    recvstream.str(recvstring);

    for (uint32_t i = layerBegin; i < layerEnd; ++i) {
      decoded_in_he[i].unsafe_load(context, recvstream);
    }
  }

  setTimePoint("Receiver::Recv back and Serialize");
}
}  // namespace rlweOkvs
