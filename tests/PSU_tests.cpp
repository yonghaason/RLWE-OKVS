#include "RPMT_tests.h"
#include "psu.h"

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
        
    PsuSender psuSender;
    PsuReceiver psuReceiver;
    psuSender.setTimer(timer_s);
    psuReceiver.setTimer(timer_r);

    psuSender.init(n, n, prng.get());
    psuReceiver.init(n, n, prng.get());

    vector<block> FX;
    vector<block> FY;
    
    timer_s.setTimePoint("start");
    timer_r.setTimePoint("start");

    for (u64 i = 0; i < 1; ++i) {
        auto p0 = psuSender.oprf(Y, FY, socket[0]);
        auto p1 = psuReceiver.oprf(X, FX, socket[1]);

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

void psu_protocol_test(const oc::CLP& cmd)
{       
    u64 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 20));
    u64 nt = cmd.getOr("nt", 1);
    
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
        
    PsuSender psuSender;
    PsuReceiver psuReceiver;
    psuSender.setTimer(timer_s);
    psuReceiver.setTimer(timer_r);

    psuSender.init(n, n, prng.get());
    psuReceiver.init(n, n, prng.get());

    vector<block> D;
    
    timer_s.setTimePoint("start");
    timer_r.setTimePoint("start");

    for (u64 i = 0; i < 1; ++i) {
        auto p0 = psuSender.run(Y, socket[0]);
        auto p1 = psuReceiver.run(X, D, socket[1]);

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

    if (D.size() != n - real) {
        cout << "Size different" << endl;
        throw RTE_LOC; 
    }

    for (size_t i = 0; i < D.size(); i++) {
        if (setX.find(D[i]) != setX.end()) {
            throw RTE_LOC;
        }
    }

    if (cmd.isSet("v")) {
        cout << endl;
        cout << timer_s << endl;
        cout << timer_r << endl;

        std::cout << "comm " << double(socket[0].bytesSent())/ 1024 / 1024 << " + "
              << double(socket[1].bytesSent())/ 1024 / 1024 << " = "
              << double(socket[0].bytesSent() + socket[1].bytesSent()) / 1024 / 1024
              << "MB" << std::endl;
    }
}