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

void sequencing_test(const oc::CLP& cmd)
{
    u64 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 20));
    double expansion_ratio = cmd.getOr("mratio", 1.16);
    u32 span_blocks = cmd.getOr("s", 20);

    Timer timer;

    PRNG prng;
    prng.SetSeed(oc::ZeroBlock);

    SspmtSender sspmtSender;
    sspmtParams ssParams;
    ssParams.initialize(n);
    sspmtSender.init(n, n, ssParams, prng.get());
    sspmtSender.rpmt_on();

    SspmtSender sspmtSender2;
    sspmtSender2.init(n, n, ssParams, prng.get());
    sspmtSender2.rpmt_on();

    uint32_t m = ceil(expansion_ratio*n);

    vector<uint32_t> start_pos(n);
    for (size_t i = 0; i < n; i++) {
        start_pos[i] = prng.get<uint32_t>() % m;
    }

    timer.setTimePoint("start");

    sspmtSender.sequencing_with_span(start_pos, span_blocks);

    timer.setTimePoint("sequencing_with_span");

    sspmtSender2.sequencing_naive(start_pos);

    timer.setTimePoint("sequencing_naive");

    std::cout << "Span:" << sspmtSender.mNumLayers << std::endl;
    std::cout << "Naive: " << sspmtSender2.mNumLayers << std::endl;

    cout << timer << endl;
}