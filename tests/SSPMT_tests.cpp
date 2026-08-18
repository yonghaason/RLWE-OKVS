#include "SSPMT_tests.h"
#include "rlwe-okvs/sspmt.h"

#include "cryptoTools/Common/Defines.h"
#include "cryptoTools/Common/CLP.h"
#include "cryptoTools/Common/TestCollection.h"
#include "cryptoTools/Common/Timer.h"
#include "cryptoTools/Crypto/PRNG.h"
#ifdef COPROTO_ENABLE_BOOST
#include <coproto/Socket/AsioSocket.h>
#endif
#include "coproto/Socket/LocalAsyncSock.h"

#include "macoro/sync_wait.h"
#include "macoro/when_all.h"
#include "macoro/thread_pool.h"

#include <iostream>
#include <vector>

using namespace std;
using namespace oc;
using namespace rlweOkvs;

void sspmt_test(const oc::CLP& cmd)
{
    u64 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 16));
    u64 nt = cmd.getOr("nt", 1);
    u64 inter = cmd.getOr("inter", n / 2);

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

    auto socket = coproto::AsioSocket::makePair();

    socket[0].setExecutor(pool0);
    socket[1].setExecutor(pool1);

    oc::Timer timer_s, timer_r;

    SspmtSender sender;
    SspmtReceiver recver;
    sender.setTimer(timer_s);
    recver.setTimer(timer_r);
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
        u64 layout = sender.getLayoutSize();
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

void sspmt_net_test(const oc::CLP& cmd)
{
    u64 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 16));
    u64 nt = cmd.getOr("nt", 1);
    u64 inter = cmd.getOr("inter", n / 2);

    string role = cmd.getOr<string>("role", "");
    if (role != "sender" && role != "recver")
        throw oc::UnitTestSkipped(
            "sspmt_net_test needs -role sender|recver (and both processes "
            "must agree on -ip <sender-addr:port>, -nn, -inter)");
    bool isSender = (role == "sender");

    string ip = cmd.getOr<string>("ip", "127.0.0.1:1212");

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
    vector<block> X(n), Y(n);
    prng.get(X.data(), X.size());
    prng.get(Y.data(), Y.size());
    for (u64 i = 0; i < inter; ++i) Y[i] = X[i];

    macoro::thread_pool pool;
    auto e = pool.make_work();
    pool.create_threads(nt);

    auto sock = coproto::asioConnect(ip, isSender);
    sock.setExecutor(pool);

    oc::Timer timer;
    oc::BitVector share;

    if (isSender)
    {
        SspmtSender sender;
        sender.setTimer(timer);
        sender.init(n, n, params, block(555, 1) );
        auto p = sender.run(Y, share, sock);
        macoro::sync_wait(std::move(p) | macoro::start_on(pool));

        std::vector<u8> buf(share.sizeBytes());
        memcpy(buf.data(), share.data(), buf.size());
        macoro::sync_wait(sock.send(coproto::copy(buf)));

        if (cmd.isSet("v"))
        {
            u64 L = sender.getNumLayers();
            u64 layout = sender.getLayoutSize();
            std::cout << "n = " << n << ", intersection = " << inter << std::endl;
            std::cout << "L (layers) = " << L
                      << ", H (slots) = " << params.heNumSlots
                      << ", L*H = " << layout << std::endl;
        }
    }
    else
    {
        SspmtReceiver recver;
        recver.setTimer(timer);
        recver.init(n, n, params, block(555, 2) );
        auto p = recver.run(X, share, sock);
        macoro::sync_wait(std::move(p) | macoro::start_on(pool));

        std::vector<u8> buf(oc::divCeil(share.size(), 8));
        macoro::sync_wait(sock.recv(buf));
        oc::BitVector senderShare(buf.data(), share.size());

        oc::BitVector flags = share;
        flags ^= senderShare;
        u64 matches = 0;
        for (u64 i = 0; i < flags.size(); ++i)
            matches += flags[i];
        if (matches != inter)
        {
            std::cout << "wrong match count: got " << matches
                      << ", expected " << inter << std::endl;
            throw RTE_LOC;
        }
        std::cout << "verified: |X n Y| = " << matches << std::endl;
    }

    macoro::sync_wait(sock.flush());

    if (cmd.isSet("v"))
    {
        std::cout << "[" << role << "]" << std::endl;
        std::cout << timer << std::endl;
        double recvByte = sock.bytesReceived();
        double sentByte = sock.bytesSent();
        std::cout << "recv " << recvByte / 1024.0 / 1024.0 << " MB, sent "
                  << sentByte / 1024.0 / 1024.0 << " MB" << std::endl;
    }
}
