#include "SEQ_tests.h"
#include "rlwe-okvs/sspmt.h"

#include "cryptoTools/Common/Defines.h"
#include "cryptoTools/Common/CLP.h"
#include "cryptoTools/Common/Timer.h"
#include "cryptoTools/Crypto/PRNG.h"

#include "seal/seal.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <set>
#include <unordered_set>

using namespace std;
using namespace oc;
using namespace seal;
using namespace rlweOkvs;

#include <stdexcept>

// Exact lower bound on the achievable layer count, via the LP dual of the
// sequencing problem: the maximum over families of block intervals J_1..J_K
// whose spanBlocks-extensions are pairwise disjoint of sum_k (max per-bin item
// count inside J_k). Every valid layer partition needs at least this many
// layers -- each layer's anchor block lies in at most one extension -- so
// equality with SspmtSender::sequenceLayers() certifies that both are optimal
// for the given instance. Test-only: the protocol never needs the bound, only
// the algorithm it certifies.
static uint64_t sequencingLowerBound(const std::vector<uint32_t> &itemBin,
                              const std::vector<uint32_t> &itemBlock,
                              uint32_t numSlots, uint32_t spanBlocks) {
  const uint32_t n = static_cast<uint32_t>(itemBin.size());
  if (n == 0) {
    return 0;
  }
  const uint32_t W = std::max<uint32_t>(spanBlocks, 1);

  uint32_t b = 0;
  for (uint32_t i = 0; i < n; ++i) {
    b = std::max(b, itemBlock[i]);
  }
  b += 1;

  std::vector<std::vector<uint32_t>> block_items(b);
  for (uint32_t i = 0; i < n; ++i) {
    block_items[itemBlock[i]].push_back(i);
  }

  // g[y] = best dual family value using intervals J = [x', y'] with y' <= y,
  // where a family is admissible if the (clipped) extensions
  // [max(0, x' - W + 1), y'] are pairwise disjoint. Taking J = [x, y] as the
  // last interval forces the previous ones to end at or before x - W.
  // Mrun[x] maintains M([x, y]) = max_bin |{items with bin, block in [x, y]}|
  // as y sweeps; cnt[x][bin] are the per-interval bin counters.
  std::vector<std::vector<uint32_t>> cnt(b,
                                         std::vector<uint32_t>(numSlots, 0));
  std::vector<uint32_t> Mrun(b, 0);
  std::vector<uint64_t> g(b, 0);

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
    uint64_t best = (y > 0) ? g[y - 1] : 0;
    for (uint32_t x = 0; x <= y; ++x) {
      if (Mrun[x] == 0) {
        continue;
      }
      const uint64_t pre = (x >= W) ? g[x - W] : 0;
      best = std::max(best, pre + Mrun[x]);
    }
    g[y] = best;
  }
  return g[b - 1];
}

// Certifies that the greedy sequencing attains the exact optimum of the
// underlying interval-covering LP. sequenceLayers() outputs a feasible
// partition and sequencingLowerBound() computes the LP-dual value, which
// lower-bounds EVERY feasible partition (weak duality); equality therefore
// proves both optimal on the instance. Part 1 sweeps randomized small
// instances, part 2 runs the production parameter sets.
void opti_sequencing_test(const oc::CLP& cmd)
{
    PRNG prng(oc::toBlock(cmd.getOr("seed", 1)));

    // --- Part 1: randomized small instances ---
    u64 trials = cmd.getOr("trials", 2000);
    for (u64 t = 0; t < trials; ++t) {
        uint32_t N = 1 + prng.get<uint32_t>() % 6;
        uint32_t b = 1 + prng.get<uint32_t>() % 12;
        uint32_t W = 1 + prng.get<uint32_t>() % 5;
        uint32_t n = prng.get<uint32_t>() % 41;

        std::vector<uint32_t> bins(n), blks(n);
        for (uint32_t i = 0; i < n; ++i) {
            bins[i] = prng.get<uint32_t>() % N;
            blks[i] = prng.get<uint32_t>() % b;
        }

        std::vector<uint32_t> itemToLayer, layerMin, layerMax;
        uint64_t L = SspmtSender::sequenceLayers(bins, blks, N, W, itemToLayer,
                                                 layerMin, layerMax);
        uint64_t LB = sequencingLowerBound(bins, blks, N, W);

        // Feasibility of the produced partition (span + one item per bin).
        auto checkFeasible = [&](uint64_t layers,
                                 const std::vector<uint32_t>& toLayer,
                                 const std::vector<uint32_t>& lmin,
                                 const std::vector<uint32_t>& lmax) {
            std::vector<std::set<uint32_t>> layerBins(layers);
            for (uint32_t i = 0; i < n; ++i) {
                uint32_t l = toLayer[i];
                if (l >= layers || layerBins[l].count(bins[i]) ||
                    blks[i] < lmin[l] || blks[i] > lmax[l] ||
                    lmax[l] - lmin[l] + 1 > W)
                    throw RTE_LOC;
                layerBins[l].insert(bins[i]);
            }
        };
        checkFeasible(L, itemToLayer, layerMin, layerMax);

        // The optimal sequencer has to be feasible too, and to hit the bound.
        std::vector<uint32_t> optToLayer, optMin, optMax;
        uint64_t Lopt = SspmtSender::sequenceLayersOptimal(
            bins, blks, N, b, W, optToLayer, optMin, optMax);
        checkFeasible(Lopt, optToLayer, optMin, optMax);
        if (Lopt != LB) {
            std::cout << "OPT-alg mismatch: N=" << N << " b=" << b << " W=" << W
                      << " n=" << n << " opt=" << Lopt << " dual=" << LB
                      << std::endl;
            throw RTE_LOC;
        }

        if (L != LB) {
            std::cout << "OPT mismatch: N=" << N << " b=" << b << " W=" << W
                      << " n=" << n << " greedy=" << L << " dual=" << LB
                      << "\nitems:";
            for (uint32_t i = 0; i < n; ++i)
                std::cout << " (" << bins[i] << "," << blks[i] << ")";
            std::cout << std::endl;
            throw RTE_LOC;
        }
    }

    // --- Part 2: production parameter sets ---
    u64 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 16));
    auto repeat = cmd.getOr("repeat", 3);

    sspmtParams params;
    params.initialize(n);
    params.bandWidth = cmd.getOr("w", params.bandWidth);
    params.bandExpansion = cmd.getOr("m_r", params.bandExpansion);
    params.span_blocks = cmd.getOr("seq_span", params.span_blocks);

    uint32_t N = params.heNumSlots;
    uint64_t m = roundUpTo(params.bandExpansion * n, N);
    uint32_t W = params.resolveSpanBlocks(n);
    // Positions are band starts, so they stop bandWidth-1 short of the end.
    uint64_t range = m - params.bandWidth + 1;
    uint32_t budget = params.resolveLayerBudget(n);

    double sumL = 0, sumMax = 0;
    for (int r = 0; r < repeat; ++r) {
        std::vector<uint32_t> bins(n), blks(n);
        for (u64 i = 0; i < n; ++i) {
            uint64_t pos = prng.get<uint64_t>() % range;
            bins[i] = (uint32_t)(pos % N);
            blks[i] = (uint32_t)(pos / N);
        }

        std::vector<uint32_t> itemToLayer, layerMin, layerMax;
        Timer timer;
        timer.setTimePoint("begin");
        uint64_t L = SspmtSender::sequenceLayers(bins, blks, N, W, itemToLayer,
                                                 layerMin, layerMax);
        timer.setTimePoint("greedy");
        std::vector<uint32_t> optToLayer, optMin, optMax;
        uint64_t Lopt = SspmtSender::sequenceLayersOptimal(
            bins, blks, N, (uint32_t)(m / N), W, optToLayer, optMin, optMax);
        timer.setTimePoint("optimal");
        uint64_t LB = sequencingLowerBound(bins, blks, N, W);
        timer.setTimePoint("dual");
        if (cmd.isSet("v")) std::cout << timer << std::endl;

        if (Lopt != LB) {
            std::cout << "OPT-alg mismatch on production-size instance: opt="
                      << Lopt << " dual=" << LB << std::endl;
            throw RTE_LOC;
        }

        uint32_t maxload = 0;
        {
            std::vector<uint32_t> cnt(N, 0);
            for (u64 i = 0; i < n; ++i) maxload = std::max(maxload, ++cnt[bins[i]]);
        }

        std::cout << "  run " << r << ": greedy L=" << L << ", optimal L="
                  << Lopt << ", dual LB=" << LB << ", max bin load=" << maxload
                  << std::endl;

        if (L != LB) {
            std::cout << "OPT mismatch on production-size instance" << std::endl;
            throw RTE_LOC;
        }
        if (L > budget) {
            std::cout << "layer budget " << budget << " exceeded by " << L
                      << std::endl;
            throw RTE_LOC;
        }
        sumL += (double)L;
        sumMax += (double)maxload;
    }

    // What the padding costs: the transmitted layer count is the public
    // budget, the realized one is what the optimum needs.
    const double avgL = sumL / repeat;
    std::cout << "n=" << n << " b=" << (m / N) << " W=" << W << " (w="
              << params.bandWidth << ", m/n=" << params.bandExpansion << ")"
              << "\n  E[L] = " << avgL
              << ", E[max bin load] = " << sumMax / repeat
              << "\n  layer budget = " << budget << ", slack = "
              << (double)budget / avgL << "x" << std::endl;
}
