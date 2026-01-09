#include "RPMT_tests.h"
#include "rlwe-okvs/sspmt.h"

#include "cryptoTools/Common/Defines.h"
#include "cryptoTools/Common/CLP.h"
#include "cryptoTools/Common/Timer.h"
#include "cryptoTools/Crypto/PRNG.h"

#include "seal/seal.h"

#include <vector>
#include <cmath>
#include <iostream>
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
        SspmtSender sspmtSender;
        sspmtParams ssParams;
        ssParams.initialize(n);
        ssParams.span_blocks = span_blocks;
        ssParams.bandExpansion = m_ratio;
        ssParams.bandWidth = w;

        sspmtSender.init(n, n, ssParams, prng.get());
        sspmtSender.rpmt_on();

        uint32_t m = ceil(ssParams.bandExpansion * n);

        vector<uint32_t> start_pos(n);
        for (size_t i = 0; i < n; i++) {
            start_pos[i] = prng.get<uint32_t>() % m;
        }
      
        sspmtSender.sequencing(start_pos);

        span_avg += sspmtSender.getNumLayers();

        auto res = max_occupied_bin_modN(start_pos, ssParams.heNumSlots);
        maxbin_avg += res.max_load;
    }
    
    timer.setTimePoint("sequencing_with_span");

    
    if (cmd.isSet("v")) {
        std::cout << "Span:" << (double) span_avg / repeat
        <<" / Maximal bin size: " << (double) maxbin_avg / repeat
        << std::endl;

        cout << timer << endl;
    }
}