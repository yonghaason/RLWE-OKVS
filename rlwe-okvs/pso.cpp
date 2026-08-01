#include "pso.h"

#include <memory>
#include <set>

#include "GMW/Gmw.h"

using namespace std;
using namespace seal;
using namespace oc;
using namespace volePSI;

namespace rlweOkvs
{
namespace {
// Bit width of the arithmetic shares produced by weightedSum*. Cardinalities
// and payload sums are reduced mod 2^32.
constexpr u64 kShareBits = 32;

// Number of GMW instances used for the (single) threshold comparison. The
// circuit is evaluated on one row of real input; the rest are padding so the
// silent OT backend stays in a comfortable size regime.
constexpr u64 kThresholdInstances = 128;
}  // namespace

Proto weightedSumSender(SilentOtExtSender &ot, PRNG &prng, Socket &chl,
                        const BitVector &shareBits,
                        const std::vector<std::vector<u32>> &weights,
                        std::vector<u32> &shares)
{
    const u64 nSlots = shareBits.size();
    const u64 k = weights.size();

    std::vector<std::array<block, 2>> rot(nSlots * k);
    co_await ot.send(rot, prng, chl);

    // The receiver holds rot[i][u]; sending delta = r0 - r1 + (1-2v)w lets it
    // reach r0 + u(1-2v)w, whose distance to the sender's r0 - vw is exactly
    // (u ^ v) * w.
    std::vector<u32> delta(nSlots * k);
    shares.assign(k, 0);
    for (u64 t = 0; t < k; ++t) {
        for (u64 s = 0; s < nSlots; ++s) {
            const u64 i = t * nSlots + s;
            const u32 r0 = (u32)rot[i][0].mData[0];
            const u32 r1 = (u32)rot[i][1].mData[0];
            const u32 w = weights[t][s];
            const u32 v = shareBits[s];
            delta[i] = r0 - r1 + (v ? (u32)(0u - w) : w);
            shares[t] += r0 - (v ? w : 0u);
        }
    }

    co_await chl.send(std::move(delta));
}

Proto weightedSumReceiver(SilentOtExtReceiver &ot, PRNG &prng, Socket &chl,
                          const BitVector &shareBits, u64 numWeights,
                          std::vector<u32> &shares)
{
    const u64 nSlots = shareBits.size();
    const u64 k = numWeights;

    BitVector choices(nSlots * k);
    for (u64 t = 0; t < k; ++t) {
        for (u64 s = 0; s < nSlots; ++s) {
            choices[t * nSlots + s] = shareBits[s];
        }
    }

    std::vector<block> rot(nSlots * k);
    co_await ot.receive(choices, rot, prng, chl);

    std::vector<u32> delta;
    co_await chl.recvResize(delta);

    shares.assign(k, 0);
    for (u64 t = 0; t < k; ++t) {
        for (u64 s = 0; s < nSlots; ++s) {
            const u64 i = t * nSlots + s;
            shares[t] += (u32)rot[i].mData[0] + (shareBits[s] ? delta[i] : 0u);
        }
    }
}

Proto PsuSender::setup(Socket &chl)
{
    if (mSetupDone) {
        co_return;
    }
    sspmtSender.init(mN, mNother, mSsParams, mPrng.get());
    sspmtSender.setTimer(getTimer());
    co_await sspmtSender.setup(chl);
    mSetupDone = true;
};

Proto PsuSender::run(const std::vector<oc::block> &Y, Socket &chl)
{
    co_await setup(chl);

    BitVector sspmt;
    co_await sspmtSender.run(Y, sspmt, chl);

    const auto &slotToItem = sspmtSender.getSlotToItem();
    const u64 nSlots = sspmt.size();

    auto comm = chl.bytesSent() + chl.bytesReceived();

    vector<array<block, 2>> rotMsgs(nSlots);
    co_await otSender.send(rotMsgs, mPrng, chl);

    // Send only the message the receiver recovers when the shares agree; if
    // they disagree it xors with the other pad and gets garbage. Empty slots
    // carry garbage too, which the receiver drops by the same test.
    vector<block> otp(nSlots);
    for (u64 s = 0; s < nSlots; s++) {
        const uint32_t item = slotToItem[s];
        const block payload = (item == UINT32_MAX) ? mPrng.get<block>() : Y[item];
        otp[s] = rotMsgs[s][sspmt[s]] ^ payload;
    }

    co_await chl.send(std::move(otp));

    comm = chl.bytesSent() + chl.bytesReceived() - comm;
    cout << "FinalOT(PSU) takes " << comm << " bytes" << endl;

    setTimePoint("Sender::Final OT PSU");
};

Proto PsuReceiver::setup(Socket &chl)
{
    if (mSetupDone) {
        co_return;
    }
    sspmtReceiver.init(mN, mNother, mSsParams, mPrng.get());
    sspmtReceiver.setTimer(getTimer());
    co_await sspmtReceiver.setup(chl);
    mSetupDone = true;
};

Proto PsuReceiver::run(const std::vector<oc::block> &X,
                       std::vector<oc::block> &D, Socket &chl)
{
    co_await setup(chl);

    BitVector sspmt;
    co_await sspmtReceiver.run(X, sspmt, chl);

    vector<block> rotMsgs(sspmt.size());
    co_await otReceiver.receive(sspmt, rotMsgs, mPrng, chl);

    vector<block> otp;
    co_await chl.recvResize(otp);

    for (size_t i = 0; i < otp.size(); i++) {
        auto tmp = otp[i] ^ rotMsgs[i];
        if (tmp.mData[1] == 0) {
            D.push_back(tmp);
        }
    }

    setTimePoint("Receiver::Final OT PSU");
};

Proto PsiCardSender::setup(Socket &chl)
{
    if (mSetupDone) {
        co_return;
    }
    sspmtSender.init(mN, mNother, mSsParams, mPrng.get());
    sspmtSender.setTimer(getTimer());
    co_await sspmtSender.setup(chl);
    mSetupDone = true;
};

Proto PsiCardSender::run(const std::vector<oc::block> &Y, Socket &chl)
{
    co_await setup(chl);

    BitVector sspmt;
    co_await sspmtSender.run(Y, sspmt, chl);

    auto comm = chl.bytesSent() + chl.bytesReceived();

    std::vector<std::vector<u32>> weights(1, std::vector<u32>(sspmt.size(), 1));
    std::vector<u32> shares;
    co_await weightedSumSender(otSender, mPrng, chl, sspmt, weights, shares);

    // Only the receiver learns the cardinality.
    co_await chl.send(shares[0]);

    comm = chl.bytesSent() + chl.bytesReceived() - comm;
    cout << "Cardinality phase takes " << comm << " bytes" << endl;

    setTimePoint("Sender::Cardinality");
};

Proto PsiCardReceiver::setup(Socket &chl)
{
    if (mSetupDone) {
        co_return;
    }
    sspmtReceiver.init(mN, mNother, mSsParams, mPrng.get());
    sspmtReceiver.setTimer(getTimer());
    co_await sspmtReceiver.setup(chl);
    mSetupDone = true;
};

Proto PsiCardReceiver::run(const std::vector<oc::block> &X,
                           oc::u64 &cardinality, Socket &chl)
{
    co_await setup(chl);

    BitVector sspmt;
    co_await sspmtReceiver.run(X, sspmt, chl);

    std::vector<u32> shares;
    co_await weightedSumReceiver(otReceiver, mPrng, chl, sspmt, 1, shares);

    u32 senderShare;
    co_await chl.recv(senderShare);

    cardinality = (u32)(shares[0] - senderShare);
    setTimePoint("Receiver::Cardinality");
};

Proto PsiCardSumSender::setup(Socket &chl)
{
    if (mSetupDone) {
        co_return;
    }
    sspmtSender.init(mN, mNother, mSsParams, mPrng.get());
    sspmtSender.setTimer(getTimer());
    co_await sspmtSender.setup(chl);
    mSetupDone = true;
};

Proto PsiCardSumSender::run(const std::vector<oc::block> &Y,
                            const std::vector<oc::u32> &payloads, Socket &chl)
{
    if (payloads.size() != Y.size()) {
        throw RTE_LOC;
    }

    co_await setup(chl);

    BitVector sspmt;
    co_await sspmtSender.run(Y, sspmt, chl);

    const auto &slotToItem = sspmtSender.getSlotToItem();
    const u64 nSlots = sspmt.size();

    auto comm = chl.bytesSent() + chl.bytesReceived();

    // One weight vector for the cardinality, one carrying the payload of the
    // item at each slot (empty slots contribute nothing).
    std::vector<std::vector<u32>> weights(2);
    weights[0].assign(nSlots, 1);
    weights[1].resize(nSlots);
    for (u64 s = 0; s < nSlots; ++s) {
        const uint32_t item = slotToItem[s];
        weights[1][s] = (item == UINT32_MAX) ? 0u : payloads[item];
    }

    std::vector<u32> shares;
    co_await weightedSumSender(otSender, mPrng, chl, sspmt, weights, shares);
    co_await chl.send(std::move(shares));

    comm = chl.bytesSent() + chl.bytesReceived() - comm;
    cout << "Card-sum phase takes " << comm << " bytes" << endl;

    setTimePoint("Sender::Card Sum");
};

Proto PsiCardSumReceiver::setup(Socket &chl)
{
    if (mSetupDone) {
        co_return;
    }
    sspmtReceiver.init(mN, mNother, mSsParams, mPrng.get());
    sspmtReceiver.setTimer(getTimer());
    co_await sspmtReceiver.setup(chl);
    mSetupDone = true;
};

Proto PsiCardSumReceiver::run(const std::vector<oc::block> &X,
                              oc::u64 &cardinality, oc::u64 &payloadSum,
                              Socket &chl)
{
    co_await setup(chl);

    BitVector sspmt;
    co_await sspmtReceiver.run(X, sspmt, chl);

    std::vector<u32> shares;
    co_await weightedSumReceiver(otReceiver, mPrng, chl, sspmt, 2, shares);

    std::vector<u32> senderShares;
    co_await chl.recvResize(senderShares);

    cardinality = (u32)(shares[0] - senderShares[0]);
    payloadSum = (u32)(shares[1] - senderShares[1]);
    setTimePoint("Receiver::Card Sum");
};

Proto PsiThresholdSender::setup(Socket &chl)
{
    if (mSetupDone) {
        co_return;
    }
    sspmtSender.init(mN, mNother, mSsParams, mPrng.get());
    sspmtSender.setTimer(getTimer());
    co_await sspmtSender.setup(chl);
    mSetupDone = true;
};

Proto PsiThresholdSender::run(const std::vector<oc::block> &Y, oc::u32 threshold,
                              Socket &chl)
{
    // "at least 0" is vacuously true and would underflow the strict form.
    if (threshold == 0) {
        throw RTE_LOC;
    }

    co_await setup(chl);

    BitVector sspmt;
    co_await sspmtSender.run(Y, sspmt, chl);

    auto comm = chl.bytesSent() + chl.bytesReceived();

    std::vector<std::vector<u32>> weights(1, std::vector<u32>(sspmt.size(), 1));
    std::vector<u32> shares;
    co_await weightedSumSender(otSender, mPrng, chl, sspmt, weights, shares);

    // The cardinality stays shared: the comparison against the threshold runs
    // inside GMW over (receiverShare - senderShare), so neither party sees the
    // count. Bundle 0 takes the receiver's share, bundle 1 the negated
    // sender's share and bundle 2 the public threshold; each party feeds its
    // own bundles and zeroes the others. The circuit tests the strict
    // "> threshold - 1", which is "at least threshold".
    auto cir = sumThresholdCircuit(kShareBits);
    Gmw gmw;
    gmw.setTimer(getTimer());
    gmw.init(kThresholdInstances, cir, 1, 1ull << 16, 0, mPrng.get());

    const u32 negShare = (u32)(0u - shares[0]);
    oc::Matrix<u32> negShareIn(kThresholdInstances, 1);
    oc::Matrix<u32> thresholdIn(kThresholdInstances, 1);
    for (u64 i = 0; i < kThresholdInstances; ++i) {
        negShareIn(i, 0) = negShare;
        thresholdIn(i, 0) = threshold - 1;
    }

    gmw.setZeroInput(0);
    gmw.setInput(1, negShareIn);
    gmw.setInput(2, thresholdIn);
    co_await gmw.run(chl);

    auto out = gmw.getOutputView(0);
    u8 outShare = out(0, 0) & 1;
    co_await chl.send(outShare);

    comm = chl.bytesSent() + chl.bytesReceived() - comm;
    cout << "Threshold phase takes " << comm << " bytes" << endl;

    setTimePoint("Sender::Threshold");
};

Proto PsiThresholdReceiver::setup(Socket &chl)
{
    if (mSetupDone) {
        co_return;
    }
    sspmtReceiver.init(mN, mNother, mSsParams, mPrng.get());
    sspmtReceiver.setTimer(getTimer());
    co_await sspmtReceiver.setup(chl);
    mSetupDone = true;
};

Proto PsiThresholdReceiver::run(const std::vector<oc::block> &X,
                                oc::u32 threshold, bool &aboveThreshold,
                                Socket &chl)
{
    if (threshold == 0) {
        throw RTE_LOC;
    }

    co_await setup(chl);

    BitVector sspmt;
    co_await sspmtReceiver.run(X, sspmt, chl);

    std::vector<u32> shares;
    co_await weightedSumReceiver(otReceiver, mPrng, chl, sspmt, 1, shares);

    auto cir = sumThresholdCircuit(kShareBits);
    Gmw gmw;
    gmw.setTimer(getTimer());
    gmw.init(kThresholdInstances, cir, 1, 1ull << 16, 1, mPrng.get());

    oc::Matrix<u32> shareIn(kThresholdInstances, 1);
    for (u64 i = 0; i < kThresholdInstances; ++i) {
        shareIn(i, 0) = shares[0];
    }

    gmw.setInput(0, shareIn);
    gmw.setZeroInput(1);
    gmw.setZeroInput(2);
    co_await gmw.run(chl);

    auto out = gmw.getOutputView(0);
    u8 senderShare;
    co_await chl.recv(senderShare);

    aboveThreshold = ((out(0, 0) ^ senderShare) & 1) != 0;
    setTimePoint("Receiver::Threshold");
};
}  // namespace rlweOkvs
