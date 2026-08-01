#include "SEQ_tests.h"
#include "RPMT_tests.h"
#include "rlwe-okvs/rpmt.h"
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

struct MaxBinLoadResult {
    uint32_t max_load;   
    uint32_t argmax_bin; 
    std::vector<uint32_t> bin_counts;
};

MaxBinLoadResult max_occupied_bin_modN(const std::vector<uint32_t>& start_pos, uint32_t N) {
    if (N == 0) {
        throw std::invalid_argument("N must be nonzero");
    }

    std::vector<uint32_t> counts(N, 0);

    for (uint32_t a : start_pos) {
        uint32_t bin = a % N;
        counts[bin]++;
    }

    uint32_t max_load = 0;
    uint32_t argmax_bin = 0;
    for (uint32_t bin = 0; bin < N; ++bin) {
        if (counts[bin] > max_load) {
            max_load = counts[bin];
            argmax_bin = bin;
        }
    }

    return {max_load, argmax_bin, std::move(counts)};
}


void sequencing_test(const oc::CLP& cmd)
{
    u64 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 20));
    u64 w = cmd.getOr("w", 109);
    double m_ratio = cmd.getOr("m_r", 1.2);
    u64 span_blocks = cmd.getOr("s", 20);
    auto repeat = cmd.getOr("repeat", 1);

    Timer timer;

    PRNG prng;
    prng.SetSeed(oc::ZeroBlock);

    

    size_t span_avg = 0;
    size_t maxbin_avg = 0;

    timer.setTimePoint("start");

    for (size_t r = 0; r < repeat; r++) {
        RpmtSender rpmtSender;
        rpmtParams ssParams;
        ssParams.initialize(n);
        ssParams.span_blocks = span_blocks;
        ssParams.bandExpansion = m_ratio;
        ssParams.bandWidth = w;

        rpmtSender.init(n, n, ssParams, prng.get());
        rpmtSender.sharedOutputOn();

        uint32_t m = ceil(ssParams.bandExpansion * n);

        vector<uint32_t> start_pos(n);
        for (size_t i = 0; i < n; i++) {
            start_pos[i] = prng.get<uint32_t>() % m;
        }
      
        rpmtSender.sequencing(start_pos);

        span_avg += rpmtSender.getNumLayers();

        auto res = max_occupied_bin_modN(start_pos, ssParams.heNumSlots);
        maxbin_avg += res.max_load;
    }
    
    timer.setTimePoint("sequencing_with_span");
    
    std::cout << "Span:" << (double) span_avg / repeat
    <<" / Maximal bin size: " << (double) maxbin_avg / repeat
    << std::endl;
    
    if (cmd.isSet("v")) {
        cout << timer << endl;
    }
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
        uint64_t L = sequenceLayers(bins, blks, N, W, itemToLayer, layerMin,
                                    layerMax);
        uint64_t LB = sequencingLowerBound(bins, blks, N, W);

        // Feasibility of the produced partition (span + one item per bin).
        std::vector<std::set<uint32_t>> layerBins(L);
        for (uint32_t i = 0; i < n; ++i) {
            uint32_t l = itemToLayer[i];
            if (l >= L || layerBins[l].count(bins[i]) ||
                blks[i] < layerMin[l] || blks[i] > layerMax[l] ||
                layerMax[l] - layerMin[l] + 1 > W)
                throw RTE_LOC;
            layerBins[l].insert(bins[i]);
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
    uint32_t N = params.heNumSlots;
    uint64_t m =
        ((uint64_t)std::ceil(params.bandExpansion * n) + N - 1) / N * N;
    uint32_t W = params.span_blocks;

    for (int r = 0; r < repeat; ++r) {
        std::vector<uint32_t> bins(n), blks(n);
        for (u64 i = 0; i < n; ++i) {
            uint64_t pos = prng.get<uint64_t>() % m;
            bins[i] = (uint32_t)(pos % N);
            blks[i] = (uint32_t)(pos / N);
        }

        std::vector<uint32_t> itemToLayer, layerMin, layerMax;
        uint64_t L = sequenceLayers(bins, blks, N, W, itemToLayer, layerMin,
                                    layerMax);
        uint64_t LB = sequencingLowerBound(bins, blks, N, W);

        uint32_t maxload = 0;
        {
            std::vector<uint32_t> cnt(N, 0);
            for (u64 i = 0; i < n; ++i) maxload = std::max(maxload, ++cnt[bins[i]]);
        }

        std::cout << "n=" << n << " b=" << (m / N) << " W=" << W
                  << " : greedy L=" << L << ", dual LB=" << LB
                  << ", max bin load=" << maxload << std::endl;

        if (L != LB) {
            std::cout << "OPT mismatch on production-size instance" << std::endl;
            throw RTE_LOC;
        }
    }
}
