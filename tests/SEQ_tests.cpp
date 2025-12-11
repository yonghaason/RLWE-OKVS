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
    u64 w = cmd.getOr("w", 109);
    double m_ratio = cmd.getOr("m_r", 1.2);
    u64 span_blocks = cmd.getOr("s", 20);

    Timer timer;

    PRNG prng;
    prng.SetSeed(oc::ZeroBlock);

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

    timer.setTimePoint("start");

    sspmtSender.sequencing(start_pos);

    timer.setTimePoint("sequencing_with_span");

    
    if (cmd.isSet("v")) {
        std::cout << "Span:" << sspmtSender.getNumLayers()
        <<" / Maximal bin size: " << sspmtSender.getMaxbinsize() 
        << std::endl;

        cout << timer << endl;
    }
}