#include "PSO_tests.h"
#include "pso.h"
#ifdef COPROTO_ENABLE_BOOST
#include <coproto/Socket/AsioSocket.h>
#endif

#include "cryptoTools/Common/Defines.h"
#include "cryptoTools/Common/CLP.h"
#include "cryptoTools/Common/Timer.h"
#include "cryptoTools/Crypto/PRNG.h"

#include "macoro/sync_wait.h"
#include "macoro/when_all.h"
#include "macoro/thread_pool.h"

#include "seal/seal.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <unordered_set>

using namespace std;
using namespace oc;
using namespace seal;
using namespace rlweOkvs;

namespace {

// Shared setup for the PSO tests: two random sets of size n whose intersection
// has a known size, the ss-PMT parameters for that size, and a pair of thread
// pools plus a socket.
struct PsoFixture
{
    u64 n, nt, inter;
    vector<block> X, Y;
    sspmtParams params;
    PRNG prng;

    macoro::thread_pool pool0, pool1;
    optional<macoro::thread_pool::work> e0, e1;
    array<coproto::AsioSocket, 2> socket;
    Timer timer_s, timer_r;
    // Bytes on the wire when the offline phase ended, so the report can split
    // setup from the online protocol.
    u64 setupSent = 0, setupRecv = 0;

    PsoFixture(const oc::CLP& cmd)
        : prng(block(9871234, 1276353))
        , socket(coproto::AsioSocket::makePair())
    {
        n = cmd.getOr("n", 1ull << cmd.getOr("nn", 16));
        nt = cmd.getOr("nt", 1);
        inter = cmd.getOr("inter", n / 2);

        // The HE parameter table only has entries for these sizes.
        u64 paramN = (n == (1ull << 16) || n == (1ull << 18) ||
                      n == (1ull << 20) || n == (1ull << 22))
                         ? n
                         : (1ull << 16);
        params.initialize(paramN);
        params.bandWidth = cmd.getOr("w", params.bandWidth);
        params.bandExpansion = cmd.getOr("m_r", params.bandExpansion);
        params.span_blocks = cmd.getOr("seq_span", params.span_blocks);
        params.layerBudget = cmd.getOr("L", params.layerBudget);

        // Items live in the low 64 bits: the PSU transfer marks garbage by a
        // nonzero high word.
        X.resize(n);
        Y.resize(n);
        for (u64 i = 0; i < n; ++i) {
            X[i] = block(0, prng.get<u64>());
            Y[i] = block(0, prng.get<u64>());
        }
        for (u64 i = 0; i < inter; ++i) Y[i] = X[i];

        e0 = pool0.make_work();
        pool0.create_threads(nt);
        e1 = pool1.make_work();
        pool1.create_threads(nt);
        socket[0].setExecutor(pool0);
        socket[1].setExecutor(pool1);
    }

    template <typename P0, typename P1>
    void runBoth(P0&& p0, P1&& p1)
    {
        auto r = macoro::sync_wait(macoro::when_all_ready(
            std::move(p0) | macoro::start_on(pool0),
            std::move(p1) | macoro::start_on(pool1)));
        std::get<0>(r).result();
        std::get<1>(r).result();
    }

    // Offline phase: the correlated randomness, which depends only on the
    // public parameters. Run separately so its cost is attributed on its own.
    template <typename S, typename R>
    void runSetup(S& sender, R& receiver)
    {
        timer_s.setTimePoint("start");
        timer_r.setTimePoint("start");
        runBoth(sender.setup(socket[0]), receiver.setup(socket[1]));
        setupSent = socket[0].bytesSent();
        setupRecv = socket[1].bytesSent();
    }

    void report(const oc::CLP& cmd)
    {
        if (!cmd.isSet("v")) return;
        auto numslots = params.heNumSlots;
        auto m = roundUpTo(params.bandExpansion * n, numslots);
        cout << "\n-------ssPMT Params-------" << endl;
        cout << "n: " << n << ", intersection: " << inter << endl;
        cout << "w: " << params.bandWidth << ", m/n: " << params.bandExpansion
             << ", seq_span: " << params.resolveSpanBlocks(n) << endl;
        cout << "layer budget: " << params.resolveLayerBudget(n)
             << " (blocks b = " << m / numslots << ")" << endl;
        cout << "--------------------------" << endl;
        cout << timer_s << endl;
        cout << timer_r << endl;
        const double MB = 1024 * 1024;
        const double offline = double(setupSent + setupRecv) / MB;
        const double total = double(socket[0].bytesSent() +
                                    socket[1].bytesSent()) / MB;
        cout << "comm: offline " << offline << " MB + online "
             << (total - offline) << " MB = " << total << " MB"
             << "   (sender->receiver "
             << double(socket[0].bytesSent()) / MB << ", receiver->sender "
             << double(socket[1].bytesSent()) / MB << ")" << endl;
    }
};

}  // namespace

void psu_test(const oc::CLP& cmd)
{
    PsoFixture f(cmd);

    PsuSender psuSender;
    PsuReceiver psuReceiver;
    psuSender.setTimer(f.timer_s);
    psuReceiver.setTimer(f.timer_r);
    psuSender.initWithParam(f.n, f.n, f.params, f.prng.get());
    psuReceiver.initWithParam(f.n, f.n, f.params, f.prng.get());

    f.runSetup(psuSender, psuReceiver);

    vector<block> D;
    f.runBoth(psuSender.run(f.Y, f.socket[0]),
              psuReceiver.run(f.X, D, f.socket[1]));

    std::unordered_set<oc::block> setX(f.X.begin(), f.X.end());
    u64 real = 0;
    for (const auto& y : f.Y)
        if (setX.count(y)) ++real;

    if (D.size() != f.n - real) {
        cout << "Size different: " << D.size() << " (proto) v.s. "
             << f.n - real << " (real)" << endl;
        throw RTE_LOC;
    }
    for (const auto& d : D)
        if (setX.count(d)) throw RTE_LOC;

    f.report(cmd);
}

void psi_card_test(const oc::CLP& cmd)
{
    PsoFixture f(cmd);

    PsiCardSender sender;
    PsiCardReceiver receiver;
    sender.setTimer(f.timer_s);
    receiver.setTimer(f.timer_r);
    sender.initWithParam(f.n, f.n, f.params, f.prng.get());
    receiver.initWithParam(f.n, f.n, f.params, f.prng.get());

    f.runSetup(sender, receiver);

    u64 cardinality = 0;
    f.runBoth(sender.run(f.Y, f.socket[0]),
              receiver.run(f.X, cardinality, f.socket[1]));

    if (cardinality != f.inter) {
        cerr << "Protocol = " << cardinality << ", real = " << f.inter << endl;
        throw RTE_LOC;
    }
    if (cmd.isSet("v")) cout << "PSI cardinality = " << cardinality << endl;

    f.report(cmd);
}

void psi_sum_test(const oc::CLP& cmd)
{
    PsoFixture f(cmd);

    vector<u32> payloads(f.n);
    for (auto& p : payloads) p = f.prng.get<u32>() % 1000;

    PsiSumSender sender;
    PsiSumReceiver receiver;
    sender.setTimer(f.timer_s);
    receiver.setTimer(f.timer_r);
    sender.initWithParam(f.n, f.n, f.params, f.prng.get());
    receiver.initWithParam(f.n, f.n, f.params, f.prng.get());

    f.runSetup(sender, receiver);

    u64 payloadSum = 0;
    f.runBoth(sender.run(f.Y, payloads, f.socket[0]),
              receiver.run(f.X, payloadSum, f.socket[1]));

    std::unordered_set<oc::block> setX(f.X.begin(), f.X.end());
    u64 expectedSum = 0;
    for (u64 i = 0; i < f.n; ++i)
        if (setX.count(f.Y[i])) expectedSum += payloads[i];

    if (payloadSum != expectedSum) {
        cerr << "Protocol = " << payloadSum << ", real = " << expectedSum
             << endl;
        throw RTE_LOC;
    }
    if (cmd.isSet("v")) cout << "PSI sum = " << payloadSum << endl;

    f.report(cmd);
}

void psi_threshold_test(const oc::CLP& cmd)
{
    PsoFixture f(cmd);

    const u32 t = cmd.getOr("t", (u32)f.inter);

    PsiThresholdSender sender;
    PsiThresholdReceiver receiver;
    sender.setTimer(f.timer_s);
    receiver.setTimer(f.timer_r);
    sender.initWithParam(f.n, f.n, f.params, f.prng.get());
    receiver.initWithParam(f.n, f.n, f.params, f.prng.get());

    f.runSetup(sender, receiver);

    bool above = false;
    f.runBoth(sender.run(f.Y, t, f.socket[0]),
              receiver.run(f.X, t, above, f.socket[1]));

    const bool expected = f.inter >= t;
    if (above != expected) {
        cerr << "threshold t = " << t << ": protocol = " << above
             << ", real = " << expected << " (|X n Y| = " << f.inter << ")"
             << endl;
        throw RTE_LOC;
    }
    if (cmd.isSet("v"))
        cout << "threshold t = " << t << " -> " << above << endl;

    f.report(cmd);
}
