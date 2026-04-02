#include "RPMT_tests.h"
#include "pso.h"
#ifdef COPROTO_ENABLE_BOOST
#include <coproto/Socket/AsioSocket.h>
#endif

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

void psu_protocol_test(const oc::CLP& cmd)
{       
    u64 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 18));
    u64 nt = cmd.getOr("nt", 1);
    
    PRNG prng;
    prng.SetSeed(oc::ZeroBlock);

    vector<u64> Xorig(n);
    vector<u64> Yorig(n);
    
    prng.get(Xorig.data(), n);
    prng.get(Yorig.data(), n);

    vector<block> X(n);
    vector<block> Y(n);
    for (size_t i = 0; i < n; i++) {
        X[i].mData[0] = Xorig[i];
        Y[i].mData[0] = Yorig[i];
    }

    u64 k = prng.get<u64>() % n; //intersection size

    std::unordered_set<size_t> selX, selY;
    selX.reserve(k * 2);
    selY.reserve(k * 2);
    while (selX.size() < k) selX.insert(prng.get<u64>() % n);
    while (selY.size() < k) selY.insert(prng.get<u64>() % n);

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

    auto socket = coproto::AsioSocket::makePair();
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

    if (cmd.isSet("v")) {
        rpmtParams parms;
        parms.initialize(n);
        cout << "\n-------Params-------" << endl;
        cout << "w: " << parms.bandWidth << endl;
        cout << "m/n: " << parms.bandExpansion << endl;
        auto numslots = parms.heNumSlots;
        auto m = roundUpTo(parms.bandExpansion * n, numslots);
        cout << "wrap: " << divCeil(parms.bandWidth * numslots, m) + 1 << endl;
        cout << "seq_span: " << parms.span_blocks<< endl;
        cout << "--------------------" << endl;
    }

    vector<block> D;
    
    timer_s.setTimePoint("start");
    timer_r.setTimePoint("start");

    auto p0 = psuSender.run(Y, socket[0]);
    auto p1 = psuReceiver.run(X, D, socket[1]);

    auto r = macoro::sync_wait(
        macoro::when_all_ready(std::move(p0) | macoro::start_on(pool0),
                            std::move(p1) | macoro::start_on(pool1)));
    std::get<0>(r).result();
    std::get<1>(r).result();
   
    u64 real = 0;
    std::unordered_set<oc::block> setX(X.begin(), X.end());
    for (const auto& y : Y)
        if (setX.find(y) != setX.end()) ++real;

    if (D.size() != n - real) {
        cout << "Size different: " 
        << D.size() << " (proto) v.s. " << n - real << " (real)"  << endl;
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

void psi_card_test(const oc::CLP& cmd)
{       
    u64 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 20));
    u64 nt = cmd.getOr("nt", 1);

    rpmtParams pmtParams;
    pmtParams.initialize(n);
    
    if (cmd.isSet("v")) {
        cout << "\n-------Params-------" << endl;
        cout << "w: " << pmtParams.bandWidth << endl;
        cout << "m/n: " << pmtParams.bandExpansion << endl;
        auto numslots = pmtParams.heNumSlots;
        auto m = roundUpTo(pmtParams.bandExpansion * n, numslots);
        cout << "wrap: " << divCeil(pmtParams.bandWidth * numslots, m) + 1 << endl;
        cout << "seq_span: " << pmtParams.span_blocks<< endl;
        cout << "--------------------" << endl;
    }
    
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

    auto socket = coproto::AsioSocket::makePair();
    // auto socket = coproto::LocalAsyncSocket::makePair();
    socket[0].setExecutor(pool0);
    socket[1].setExecutor(pool1);
    
    oc::Timer timer_s;
    oc::Timer timer_r;
        
    RpmtSender rpmtSender;
    RpmtReceiver rpmtReceiver;
    rpmtSender.setTimer(timer_s);
    rpmtReceiver.setTimer(timer_r);

    rpmtSender.init(n, n, pmtParams, prng.get());
    rpmtReceiver.init(n, n, pmtParams, prng.get());

    oc::BitVector results;
    
    timer_s.setTimePoint("start");
    timer_r.setTimePoint("start");

    auto p0 = rpmtSender.run(Y, socket[0]);
    auto p1 = rpmtReceiver.run(X, results, socket[1]);

    auto r = macoro::sync_wait(
        macoro::when_all_ready(std::move(p0) | macoro::start_on(pool0),
                            std::move(p1) | macoro::start_on(pool1)));
    std::get<0>(r).result();
    std::get<1>(r).result();
   
    u64 real = 0;
    std::unordered_set<oc::block> setX(X.begin(), X.end());
    for (const auto& y : Y)
        if (setX.find(y) != setX.end()) ++real;

    u64 psi_cardinality = results.hammingWeight();

    if (psi_cardinality != real) {
        cerr << "Protocol = " << psi_cardinality
             << " / Expected = " << real << endl;
        throw RTE_LOC; 
    }

    if (cmd.isSet("v")) {
        cout << "PSI cardinality = " << psi_cardinality << endl;
        cout << "Expected        = " << real << endl;
        cout << endl;
        cout << timer_s << endl;
        cout << timer_r << endl;

        std::cout << "comm " << double(socket[0].bytesSent())/ 1024 / 1024 << " + "
              << double(socket[1].bytesSent())/ 1024 / 1024 << " = "
              << double(socket[0].bytesSent() + socket[1].bytesSent()) / 1024 / 1024
              << "MB" << std::endl;
    }
}

void psi_card_sum_32_test(const oc::CLP& cmd)
{
    u64 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 18));
    u64 nt = cmd.getOr("nt", 1);

    PRNG prng;
    prng.SetSeed(oc::ZeroBlock);
    PRNG payloadPrng;
    payloadPrng.SetSeed(oc::OneBlock);

    vector<u64> Xorig(n);
    vector<u64> Yorig(n);
    vector<u32> payloads(n);

    prng.get(Xorig.data(), n);
    prng.get(Yorig.data(), n);

    vector<block> X(n);
    vector<block> Y(n);
    for (size_t i = 0; i < n; i++) {
        X[i].mData[0] = Xorig[i];
        Y[i].mData[0] = Yorig[i];
    }

    u64 k = prng.get<u64>() % n;

    std::unordered_set<size_t> selX, selY;
    selX.reserve(k * 2);
    selY.reserve(k * 2);
    while (selX.size() < k) selX.insert(prng.get<u64>() % n);
    while (selY.size() < k) selY.insert(prng.get<u64>() % n);

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

    auto socket = coproto::AsioSocket::makePair();
    socket[0].setExecutor(pool0);
    socket[1].setExecutor(pool1);

    oc::Timer timer_s;
    oc::Timer timer_r;

    PsiCardSumSender psuSender;
    PsiCardSumReceiver psuReceiver;
    psuSender.setTimer(timer_s);
    psuReceiver.setTimer(timer_r);

    psuSender.init(n, n, prng.get());
    psuReceiver.init(n, n, prng.get());

    for (auto& payload : payloads) {
        payload = payloadPrng.get<u32>() & ((1ull << 20) - 1);
    }

    oc::u64 psiCardSum = 0;

    timer_s.setTimePoint("start");
    timer_r.setTimePoint("start");

    auto p0 = psuSender.run(Y, payloads, socket[0]);
    auto p1 = psuReceiver.run(X, psiCardSum, socket[1]);

    auto r = macoro::sync_wait(
        macoro::when_all_ready(std::move(p0) | macoro::start_on(pool0),
                            std::move(p1) | macoro::start_on(pool1)));
    std::get<0>(r).result();
    std::get<1>(r).result();

    u32 expected = 0;
    std::unordered_set<oc::block> setX(X.begin(), X.end());
    for (size_t i = 0; i < Y.size(); i++) {
        if (setX.find(Y[i]) != setX.end()) {
            expected += payloads[i];
        }
    }

    if (psiCardSum != expected) {
        cerr << "Protocol = " << psiCardSum
             << " / Expected = " << expected << endl;
        throw RTE_LOC;
    }

    if (cmd.isSet("v")) {
        cout << "PSI card sum = " << psiCardSum << endl;
        cout << "Expected     = " << expected << endl;
        cout << endl;
        cout << timer_s << endl;
        cout << timer_r << endl;

        std::cout << "comm " << double(socket[0].bytesSent()) / 1024 / 1024 << " + "
              << double(socket[1].bytesSent()) / 1024 / 1024 << " = "
              << double(socket[0].bytesSent() + socket[1].bytesSent()) / 1024 / 1024
              << "MB" << std::endl;
    }
}
