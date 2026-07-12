#include "SSPMT_tests.h"
#include "rlwe-okvs/sspmt.h"

#include "cryptoTools/Common/Defines.h"
#include "cryptoTools/Common/CLP.h"
#include "cryptoTools/Common/Timer.h"
#include "cryptoTools/Crypto/PRNG.h"
#include "coproto/Socket/LocalAsyncSock.h"

#include "macoro/sync_wait.h"
#include "macoro/when_all.h"
#include "macoro/thread_pool.h"

#include <iostream>
#include <vector>

using namespace std;
using namespace oc;
using namespace rlweOkvs;

// End-to-end correctness of the full-layout RLWE-OKVS ss-PMT.
//
// Roles: the sender holds Y and drives the sequenced homomorphic decode; the
// receiver holds X and encodes the indicator OKVS. The output is a per-slot
// XOR share over the WHOLE L x H layout (no occupancy is revealed). A slot
// reconstructs to 1 iff the Y-item sitting there is in X, and empty slots
// reconstruct to 0, so the reconstructed bits sum to exactly |X n Y| -- which
// is the ground truth we check, without needing the sender's internal layout.
void sspmt_fulllayout_protocol_test(const oc::CLP& cmd)
{
    u64 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 16));
    u64 nt = cmd.getOr("nt", 1);
    u64 inter = cmd.getOr("inter", n / 2);

    // The HE parameter table only has entries for these sizes.
    u64 paramN = (n == (1ull << 16) || n == (1ull << 18) ||
                  n == (1ull << 20) || n == (1ull << 22))
                     ? n
                     : (1ull << 16);
    sspmtParams params;
    params.initialize(paramN);
    params.bandWidth = cmd.getOr("w", params.bandWidth);
    params.bandExpansion = cmd.getOr("m_r", params.bandExpansion);
    params.span_blocks = cmd.getOr("seq_span", params.span_blocks);

    PRNG prng(block(9871234, 1276353));

    // X, Y distinct random, agreeing on exactly `inter` items.
    vector<block> X(n), Y(n);
    prng.get(X.data(), X.size());
    prng.get(Y.data(), Y.size());
    for (u64 i = 0; i < inter; ++i) Y[i] = X[i];

    macoro::thread_pool pool0;
    auto e0 = pool0.make_work();
    pool0.create_threads(nt);
    macoro::thread_pool pool1;
    auto e1 = pool1.make_work();
    pool1.create_threads(nt);

    auto socket = coproto::LocalAsyncSocket::makePair();
    socket[0].setExecutor(pool0);
    socket[1].setExecutor(pool1);

    oc::Timer timer_s, timer_r;

    SspmtSender sender;
    SspmtReceiver recver;
    sender.setTimer(timer_s);
    recver.setTimer(timer_r);
    sender.fullLayoutOn();
    recver.fullLayoutOn();
    sender.init(n, n, params, prng.get());
    recver.init(n, n, params, prng.get());

    oc::BitVector ss, rs;

    auto p0 = sender.run(Y, ss, socket[0]);
    auto p1 = recver.run(X, rs, socket[1]);
    auto r = macoro::sync_wait(
        macoro::when_all_ready(std::move(p0) | macoro::start_on(pool0),
                               std::move(p1) | macoro::start_on(pool1)));
    std::get<0>(r).result();
    std::get<1>(r).result();

    if (ss.size() != rs.size())
        throw RTE_LOC;

    // Reconstruct: every slot's bit is the XOR of the two shares.
    oc::BitVector flags = ss;
    flags ^= rs;

    u64 matches = 0;
    for (u64 i = 0; i < flags.size(); ++i)
        matches += flags[i];

    if (matches != inter)
    {
        std::cout << "wrong match count: got " << matches
                  << ", expected " << inter << std::endl;
        throw RTE_LOC;
    }

    if (cmd.isSet("v"))
    {
        u64 L = sender.getNumLayers();
        u64 layout = sender.getLayoutSize();  // L * H
        std::cout << "n = " << n << ", intersection = " << inter << std::endl;
        std::cout << "L (layers) = " << L
                  << ", H (slots) = " << params.heNumSlots
                  << ", L*H = " << layout << std::endl;
        std::cout << "equality instances / n = "
                  << static_cast<double>(layout) / n
                  << "  (this is the L*H / n_y blow-up)" << std::endl;
        std::cout << timer_s << std::endl;
        std::cout << timer_r << std::endl;
        double recvByte = socket[0].bytesReceived();
        double sentByte = socket[0].bytesSent();
        std::cout << recvByte / 1024.0 / 1024.0 << " + "
                  << sentByte / 1024.0 / 1024.0 << " = "
                  << (recvByte + sentByte) / 1024.0 / 1024.0 << " MB "
                  << std::endl;
    }
}
