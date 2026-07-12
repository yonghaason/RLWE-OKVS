#include "CPSI_tests.h"
#include "band_okvs/cpsi.h"

#include "cryptoTools/Common/Defines.h"
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

// End-to-end correctness of the CPSI-based ss-PMT:
// the receiver holds X, the sender holds Y, |X intersect Y| = inter,
// and for every receiver item the XOR of the two flag shares at its
// cuckoo bin must equal its membership bit. Unmapped (empty) bins must
// share a 0.
void cpsi_sspmt_protocol_test(const oc::CLP& cmd)
{
    u64 nr = cmd.getOr("n", 1ull << cmd.getOr("nn", 12));
    u64 ns = cmd.getOr("ns", nr);
    u64 nt = cmd.getOr("nt", 1);
    u64 inter = cmd.getOr("inter", min(nr, ns) / 2);

    PRNG prng(block(4253465, 3434565));

    // X, Y random and distinct, agreeing exactly on the first `inter` items.
    vector<block> X(nr), Y(ns);
    prng.get(X.data(), X.size());
    prng.get(Y.data(), Y.size());
    copy(X.begin(), X.begin() + inter, Y.begin());

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

    CpsiSender sender;
    CpsiReceiver recver;
    sender.setTimer(timer_s);
    recver.setTimer(timer_r);
    sender.init(ns, nr, prng.get());
    recver.init(nr, ns, prng.get());

    CpsiSharing ss, rs;

    auto p0 = sender.run(Y, ss, socket[0]);
    auto p1 = recver.run(X, rs, socket[1]);
    auto r = macoro::sync_wait(
        macoro::when_all_ready(std::move(p0) | macoro::start_on(pool0),
                               std::move(p1) | macoro::start_on(pool1)));
    std::get<0>(r).result();
    std::get<1>(r).result();

    if (ss.mNumBins != rs.mNumBins)
        throw RTE_LOC;

    // Reconstruct the flags and check them against the ground truth.
    oc::BitVector flags = ss.mFlagShares;
    flags ^= rs.mFlagShares;

    vector<bool> mapped(rs.mNumBins, false);
    for (u64 b = 0; b < nr; ++b)
    {
        auto bin = rs.mMapping[b];
        if (bin >= rs.mNumBins || mapped[bin])
            throw RTE_LOC;
        mapped[bin] = true;

        bool expected = (b < inter);
        if (flags[bin] != expected)
        {
            std::cout << "wrong flag for X[" << b << "] (bin " << bin
                      << "): got " << flags[bin]
                      << ", expected " << expected << std::endl;
            throw RTE_LOC;
        }
    }

    // Empty bins carry no item, so their shares must reconstruct to 0.
    for (u64 i = 0; i < rs.mNumBins; ++i)
    {
        if (!mapped[i] && flags[i])
        {
            std::cout << "empty bin " << i << " reconstructed to 1" << std::endl;
            throw RTE_LOC;
        }
    }

    if (cmd.isSet("v"))
    {
        std::cout << "n_r = " << nr << ", n_s = " << ns
                  << ", numBins = " << rs.mNumBins << std::endl;
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
