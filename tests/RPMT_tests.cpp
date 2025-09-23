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

void rpmt_protocol_test(const oc::CLP& cmd)
{       
    u64 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 20));
    u64 nt = cmd.getOr("nt", 1);

    // u64 logp = cmd.getOr("logp", 60);
    // u64 numSlots = cmd.getOr("t_s", 1 << 13);
    // u64 width = cmd.getOr("w", 134);
    // double expansion_ratio = cmd.getOr("mratio", 1.16);

    sspmtParams ssParams;
    // ssParams.bandExpansion = expansion_ratio;
    // ssParams.bandWidth = width;
    // ssParams.hePlainModulusBits = logp;
    // ssParams.heNumSlots = numSlots;
    
    PRNG prng;
    prng.SetSeed(oc::ZeroBlock);
    vector<block> X(n);
    vector<block> Y(n);

    prng.get(X.data(), static_cast<u64>(n));
    prng.get(Y.data(), static_cast<u64>(n));

    u64 k = 1 + (prng.get<u64>() % static_cast<u64>(n)); //intersection

    std::unordered_set<size_t> selX, selY;
    selX.reserve(k * 2);
    selY.reserve(k * 2);
    while (selX.size() < k) selX.insert(static_cast<size_t>(prng.get<u64>() % n));
    while (selY.size() < k) selY.insert(static_cast<size_t>(prng.get<u64>() % n));

    auto itX = selX.begin();
    auto itY = selY.begin();
    for (; itX != selX.end() && itY != selY.end(); ++itX, ++itY)
        Y[*itY] = X[*itX];

    macoro::thread_pool pool0;
    auto e0 = pool0.make_work();
    pool0.create_threads(nt);
    macoro::thread_pool pool1;
    auto e1 = pool1.make_work();
    pool1.create_threads(nt);

    // auto socket = coproto::AsioSocket::makePair();
    auto socket = coproto::LocalAsyncSocket::makePair();
    socket[0].setExecutor(pool0);
    socket[1].setExecutor(pool1);
    
    oc::Timer timer_s;
    oc::Timer timer_r;
        
    SspmtSender sspmtSender;
    SspmtReceiver sspmtReceiver;
    sspmtSender.setTimer(timer_s);
    sspmtReceiver.setTimer(timer_r);

    sspmtSender.init(n, n, ssParams, prng.get());
    sspmtSender.rpmt_on();
    sspmtReceiver.init(n, n, ssParams, prng.get());
    sspmtReceiver.rpmt_on();

    oc::BitVector results;
    
    timer_s.setTimePoint("start");
    timer_r.setTimePoint("start");

    for (u64 i = 0; i < 1; ++i) {
        auto p0 = sspmtSender.run(Y, socket[0]);
        auto p1 = sspmtReceiver.run(X, results, socket[1]);

        auto r = macoro::sync_wait(
            macoro::when_all_ready(std::move(p0) | macoro::start_on(pool0),
                                std::move(p1) | macoro::start_on(pool1)));
        std::get<0>(r).result();
        std::get<1>(r).result();
    }
   
    u64 real = 0;
    std::unordered_set<oc::block> setX(X.begin(), X.end());
    for (const auto& y : Y)
        if (setX.find(y) != setX.end()) ++real;

    if (results.hammingWeight() != real) {
        cerr << "Protocol = " << results.hammingWeight()
             << " / Expected = " << real << endl;
        throw RTE_LOC; 
    }

    if (cmd.isSet("v")) {
        cout << endl;
        cout << timer_s << endl;
        cout << timer_r << endl;
    } 
}