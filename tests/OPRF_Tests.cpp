#include "RPMT_tests.h"
#include "band_okvs/oprf.h"

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

void oprf_protocol_test(const oc::CLP& cmd)
{       
    u64 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 20));
    u64 nt = cmd.getOr("nt", 1);
    
    PRNG prng;
    prng.SetSeed(oc::ZeroBlock);
    vector<block> X(n);
    prng.get(X.data(), n);
    vector<block> Y = X;
    prng.get(Y.data(), 10);

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
        
    OprfSender oprfSender;
    OprfReceiver oprfReceiver;
    oprfSender.setTimer(timer_s);
    oprfReceiver.setTimer(timer_r);

    oprfSender.init(n, n, prng.get());
    oprfReceiver.init(n, n, prng.get());

    vector<block> FX;
    vector<block> FY;
    
    timer_s.setTimePoint("start");
    timer_r.setTimePoint("start");

    for (u64 i = 0; i < 1; ++i) {
        auto p0 = oprfSender.run(Y, FY, socket[0]);
        auto p1 = oprfReceiver.run(X, FX, socket[1]);

        auto r = macoro::sync_wait(
            macoro::when_all_ready(std::move(p0) | macoro::start_on(pool0),
                                std::move(p1) | macoro::start_on(pool1)));
        std::get<0>(r).result();
        std::get<1>(r).result();
    }
   
    for (size_t i = 0; i < 10; i++) {
        if (FX[i] == FY[i]) throw RTE_LOC;
    }
    
    for (size_t i = 10; i < n; i++) {
        if (FX[i] != FY[i]) throw RTE_LOC;
    }

    if (cmd.isSet("v")) {
        cout << endl;
        cout << timer_s << endl;
        cout << timer_r << endl;
    } 
}