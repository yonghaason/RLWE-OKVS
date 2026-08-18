#include "GMW_tests.h"
#include "GMW/Gmw.h"
#include "GMW/SilentTripleGen.h"
#include "cryptoTools/Network/IOService.h"
#include "cryptoTools/Network/Session.h"
#include "cryptoTools/Circuit/BetaLibrary.h"
#include "coproto/Socket/LocalAsyncSock.h"
#ifdef COPROTO_ENABLE_BOOST
#include <coproto/Socket/AsioSocket.h>
#endif
#include "cryptoTools/Common/Timer.h"
#include <numeric>

#include "macoro/task.h"
#include "macoro/sync_wait.h"
#include "macoro/when_all.h"

using coproto::LocalAsyncSocket;

using namespace volePSI;
using PRNG = oc::PRNG;

void Gmw_threshold_test(const oc::CLP& cmd)
{
    const u64 bits = 32;
    auto cir = sumThresholdCircuit(bits);

    oc::PRNG prng(oc::toBlock(cmd.getOr<u64>("s", 0)));

    auto pack = [&](u32 v) {
        oc::BitVector bv(bits);
        memcpy(bv.data(), &v, sizeof(v));
        return bv;
    };

    for (u64 i = 0; i < 200; ++i)
    {
        u32 card = prng.get<u32>() % 1000;
        u32 t = prng.get<u32>() % 1000;
        u32 share1 = prng.get<u32>();
        u32 share0 = card - share1;

        std::array<oc::BitVector, 3> in{pack(share1), pack(share0), pack(t)};
        std::array<oc::BitVector, 1> out{oc::BitVector(1)};
        cir.evaluate({in.data(), 3}, {out.data(), 1}, false);

        if (out[0][0] != (card > t))
        {
            std::cout << "plain: card " << card << " t " << t << " -> "
                      << (int)out[0][0] << std::endl;
            throw RTE_LOC;
        }
    }

    const u64 nInst = 128;
    const u32 card = 37, t = 36;
    u32 share1 = prng.get<u32>();
    u32 share0 = card - share1;

    macoro::thread_pool pool0, pool1;
    auto e0 = pool0.make_work();
    pool0.create_threads(1);
    auto e1 = pool1.make_work();
    pool1.create_threads(1);
    auto socket = coproto::LocalAsyncSocket::makePair();
    socket[0].setExecutor(pool0);
    socket[1].setExecutor(pool1);

    Gmw gmw0, gmw1;
    gmw0.init(nInst, cir, 1, 1ull << 16, 0, oc::ZeroBlock);
    gmw1.init(nInst, cir, 1, 1ull << 16, 1, oc::OneBlock);

    oc::Matrix<u32> in0(nInst, 1), inT(nInst, 1), in1(nInst, 1);
    for (u64 i = 0; i < nInst; ++i)
    {
        in0(i, 0) = share0;
        inT(i, 0) = t;
        in1(i, 0) = share1;
    }

    gmw0.setZeroInput(0);
    gmw0.setInput(1, in0);
    gmw0.setInput(2, inT);
    gmw1.setInput(0, in1);
    gmw1.setZeroInput(1);
    gmw1.setZeroInput(2);

    auto p0 = gmw0.run(socket[0]);
    auto p1 = gmw1.run(socket[1]);
    auto r = macoro::sync_wait(
        macoro::when_all_ready(std::move(p0), std::move(p1)));
    std::get<0>(r).result();
    std::get<1>(r).result();

    auto o0 = gmw0.getOutputView(0);
    auto o1 = gmw1.getOutputView(0);
    const bool got = ((o0(0, 0) ^ o1(0, 0)) & 1) != 0;
    if (got != (card > t))
    {
        std::cout << "gmw: card " << card << " t " << t << " -> " << got
                  << std::endl;
        throw RTE_LOC;
    }
}

void Gmw_iszero_test(const oc::CLP& cmd)
{
    block seed = oc::toBlock(cmd.getOr<u64>("s", 0));
    oc::PRNG prng(seed);

    u64 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 10));
    u64 bs = cmd.getOr("b", 1ull << cmd.getOr("bb", 10));
    u64 keyBitLength = cmd.getOr("kbl", 64);
    u64 nt = cmd.getOr("nt", 1);
    u64 keyByteLength = (keyBitLength + 7) / 8;

    macoro::thread_pool pool0;
    auto e0 = pool0.make_work();
    pool0.create_threads(nt);
    macoro::thread_pool pool1;
    auto e1 = pool1.make_work();
    pool1.create_threads(nt);

    auto socket = coproto::LocalAsyncSocket::makePair();
    socket[0].setExecutor(pool0);
    socket[1].setExecutor(pool1);

    Matrix<u8> in0, in1;
    in0.resize(n, keyByteLength, oc::AllocType::Uninitialized);
    in1.resize(n, keyByteLength, oc::AllocType::Uninitialized);

    prng.get<u8>(in0);
    prng.get<u8>(in1);
    auto diff = prng.get<u64>() % n;

    for (size_t i = 0; i < diff; i++) {
        memcpy(&in0(i, 0), &in1(i, 0), keyByteLength);
    }

    oc::Timer timer1, timer2;
    timer1.setTimePoint("start");
    timer2.setTimePoint("start");

    Gmw gmw0, gmw1;
    gmw0.setTimer(timer1);
    gmw1.setTimer(timer2);
    auto cir = isZeroCircuit(keyBitLength);
    gmw0.init(n, cir, 1, bs, 0, oc::ZeroBlock);
    gmw1.init(n, cir, 1, bs, 1, oc::OneBlock);

    gmw0.setInput(0, in0);
    gmw1.setInput(0, in1);

    auto p0 = gmw0.run(socket[0]);
    auto p1 = gmw1.run(socket[1]);
    auto r = macoro::sync_wait(macoro::when_all_ready(std::move(p0), std::move(p1)));
    std::get<0>(r).result();
    std::get<1>(r).result();

    auto sout0 = gmw0.getOutputView(0);
    auto sout1 = gmw1.getOutputView(0);

    oc::BitVector sout0bv(n);
    std::copy(sout0.begin(), sout0.end(), sout0bv.data());
    oc::BitVector sout1bv(n);
    std::copy(sout1.begin(), sout1.end(), sout1bv.data());

    timer1.setTimePoint("end");
    timer2.setTimePoint("end");

    for (u64 i = 0; i < diff; ++i)
    {
        auto act = sout0bv[i] ^ sout1bv[i];
        if (!act)
        {
            std::cout << "i   " << i << std::endl;
            throw RTE_LOC;
        }
    }
    for (u64 i = diff; i < n; ++i)
    {
        auto act = sout0bv[i] ^ sout1bv[i];
        if (act)
        {
            std::cout << "i   " << i << std::endl;
            throw RTE_LOC;
        }
    }

    if (cmd.isSet("v")) {
        std::cout << timer1 << std::endl;
        std::cout << timer2 << std::endl;
        double recvByte = socket[0].bytesReceived();
        double sentByte = socket[0].bytesSent();
        std::cout
        << recvByte/1024.0/1024.0 << " + "  << sentByte/1024.0/1024.0
        << " = " << (recvByte + sentByte)/1024.0/1024.0 << " MB "
        << std::endl;
    }
}
