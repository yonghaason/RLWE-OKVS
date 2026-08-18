#include "pso.h"

#include <memory>
#include <numeric>
#include <set>

#include "GMW/Gmw.h"

using namespace std;
using namespace seal;
using namespace oc;
using namespace volePSI;

namespace rlweOkvs
{
namespace {

constexpr u64 kShareBits = 32;

constexpr u64 kThresholdInstances = 128;
}

Proto genRotSender(SilentOtExtSender &ot, PRNG &prng, Socket &chl, u64 count,
                   RotStore &store)
{
    store.mSend.resize(count);
    co_await ot.silentSend(store.mSend, prng, chl);
}

Proto genRotReceiver(SilentOtExtReceiver &ot, PRNG &prng, Socket &chl,
                     u64 count, RotStore &store)
{
    store.mRecv.resize(count);
    store.mChoices.resize(count);
    co_await ot.silentReceive(store.mChoices, store.mRecv, prng, chl,
                              oc::OTType::Random);
}

Proto weightedSumSender(const RotStore &rot, Socket &chl,
                        const BitVector &shareBits,
                        const std::vector<std::vector<u32>> &weights,
                        std::vector<u32> &shares)
{
    const u64 nSlots = shareBits.size();
    const u64 k = weights.size();

    BitVector flip(nSlots * k);
    co_await chl.recv(flip);

    std::vector<u32> delta(nSlots * k);
    shares.assign(k, 0);
    for (u64 t = 0; t < k; ++t) {
        for (u64 s = 0; s < nSlots; ++s) {
            const u64 i = t * nSlots + s;
            const bool swap = flip[i];
            const u32 r0 = (u32)rot.mSend[i][swap ? 1 : 0].mData[0];
            const u32 r1 = (u32)rot.mSend[i][swap ? 0 : 1].mData[0];
            const u32 w = weights[t][s];
            const u32 v = shareBits[s];
            delta[i] = r0 - r1 + (v ? (u32)(0u - w) : w);
            shares[t] += r0 - (v ? w : 0u);
        }
    }

    co_await chl.send(std::move(delta));
}

Proto weightedSumReceiver(const RotStore &rot, Socket &chl,
                          const BitVector &shareBits, u64 numWeights,
                          std::vector<u32> &shares)
{
    const u64 nSlots = shareBits.size();
    const u64 k = numWeights;

    BitVector flip(nSlots * k);
    for (u64 t = 0; t < k; ++t) {
        for (u64 s = 0; s < nSlots; ++s) {
            flip[t * nSlots + s] = rot.mChoices[t * nSlots + s] ^ shareBits[s];
        }
    }
    co_await chl.send(flip);

    std::vector<u32> delta;
    co_await chl.recvResize(delta);

    shares.assign(k, 0);
    for (u64 t = 0; t < k; ++t) {
        for (u64 s = 0; s < nSlots; ++s) {
            const u64 i = t * nSlots + s;
            shares[t] += (u32)rot.mRecv[i].mData[0] + (shareBits[s] ? delta[i] : 0u);
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

    const u64 nSlots = sspmtSender.getLayoutSize();
    psSender.setTimer(getTimer());
    psSender.init(nSlots);
    co_await psSender.correlate(otSwitchReceiver, mPrng, chl);

    co_await genRotSender(otSender, mPrng, chl, mN, mRot);
    setTimePoint("PsuSender::Setup (transfer OTs)");

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

    std::vector<int> itemToSlot(mN, -1);
    for (u64 s = 0; s < nSlots; ++s) {
        if (slotToItem[s] != UINT32_MAX) itemToSlot[slotToItem[s]] = (int)s;
    }

    BitVector psShare;
    co_await psSender.apply(psShare, chl);

    const auto &invPerm = psSender.getInvPermRef();
    std::vector<u32> posToItem(nSlots, UINT32_MAX);
    BitVector realPos(nSlots);
    for (u64 j = 0; j < mN; ++j) {
        const u32 p = (u32)invPerm[itemToSlot[j]];
        posToItem[p] = (u32)j;
        realPos[p] = 1;
    }
    co_await chl.send(realPos);

    std::vector<u32> order(mN);
    BitVector bits(mN);
    {
        u64 k = 0;
        for (u64 p = 0; p < nSlots; ++p) {
            const u32 item = posToItem[p];
            if (item == UINT32_MAX) continue;
            order[k] = item;
            bits[k] = psShare[p] ^ sspmt[itemToSlot[item]];
            ++k;
        }
    }

    BitVector flip(mN);
    co_await chl.recv(flip);

    vector<block> otp(mN);
    for (u64 j = 0; j < mN; j++) {
        const bool swap = flip[j];
        otp[j] = mRot.mSend[j][(bits[j] ^ swap) ? 1 : 0] ^ Y[order[j]];
    }

    co_await chl.send(std::move(otp));

    comm = chl.bytesSent() + chl.bytesReceived() - comm;
    cout << "PermuteShare + FinalOT(PSU) takes " << comm << " bytes" << endl;

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

    const u64 nSlots = sspmtReceiver.getLayoutSize();
    psReceiver.setTimer(getTimer());
    psReceiver.init(nSlots);
    co_await psReceiver.correlate(otSwitchSender, mPrng, chl);

    co_await genRotReceiver(otReceiver, mPrng, chl, mNother, mRot);
    setTimePoint("PsuReceiver::Setup (transfer OTs)");

    mSetupDone = true;
};

Proto PsuReceiver::run(const std::vector<oc::block> &X,
                       std::vector<oc::block> &D, Socket &chl)
{
    co_await setup(chl);

    BitVector sspmt;
    co_await sspmtReceiver.run(X, sspmt, chl);

    BitVector psShare;
    co_await psReceiver.apply(sspmt, psShare, chl);

    BitVector realPos(psShare.size());
    co_await chl.recv(realPos);

    BitVector bits(mNother);
    {
        u64 k = 0;
        for (u64 p = 0; p < realPos.size(); ++p) {
            if (realPos[p]) bits[k++] = psShare[p];
        }
    }

    BitVector flip(mNother);
    for (u64 j = 0; j < mNother; ++j) flip[j] = mRot.mChoices[j] ^ bits[j];
    co_await chl.send(flip);

    vector<block> otp;
    co_await chl.recvResize(otp);

    for (size_t i = 0; i < otp.size(); i++) {
        auto tmp = otp[i] ^ mRot.mRecv[i];
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

    const u64 nSlots = sspmtSender.getLayoutSize();
    co_await genRotSender(otSender, mPrng, chl, nSlots, mRot);
    setTimePoint("PsiCardSender::Setup (transfer OTs)");

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
    co_await weightedSumSender(mRot, chl, sspmt, weights, shares);

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

    const u64 nSlots = sspmtReceiver.getLayoutSize();
    co_await genRotReceiver(otReceiver, mPrng, chl, nSlots, mRot);
    setTimePoint("PsiCardReceiver::Setup (transfer OTs)");

    mSetupDone = true;
};

Proto PsiCardReceiver::run(const std::vector<oc::block> &X,
                           oc::u64 &cardinality, Socket &chl)
{
    co_await setup(chl);

    BitVector sspmt;
    co_await sspmtReceiver.run(X, sspmt, chl);

    std::vector<u32> shares;
    co_await weightedSumReceiver(mRot, chl, sspmt, 1, shares);

    u32 senderShare;
    co_await chl.recv(senderShare);

    cardinality = (u32)(shares[0] - senderShare);
    setTimePoint("Receiver::Cardinality");
};

Proto PsiSumSender::setup(Socket &chl)
{
    if (mSetupDone) {
        co_return;
    }
    sspmtSender.init(mN, mNother, mSsParams, mPrng.get());
    sspmtSender.setTimer(getTimer());
    co_await sspmtSender.setup(chl);

    const u64 nSlots = sspmtSender.getLayoutSize();
    co_await genRotSender(otSender, mPrng, chl, nSlots, mRot);
    setTimePoint("PsiSumSender::Setup (transfer OTs)");

    mSetupDone = true;
};

Proto PsiSumSender::run(const std::vector<oc::block> &Y,
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

    std::vector<std::vector<u32>> weights(1);
    weights[0].resize(nSlots);
    for (u64 s = 0; s < nSlots; ++s) {
        const uint32_t item = slotToItem[s];
        weights[0][s] = (item == UINT32_MAX) ? 0u : payloads[item];
    }

    std::vector<u32> shares;
    co_await weightedSumSender(mRot, chl, sspmt, weights, shares);
    co_await chl.send(std::move(shares));

    comm = chl.bytesSent() + chl.bytesReceived() - comm;
    cout << "Sum phase takes " << comm << " bytes" << endl;

    setTimePoint("Sender::Sum");
};

Proto PsiSumReceiver::setup(Socket &chl)
{
    if (mSetupDone) {
        co_return;
    }
    sspmtReceiver.init(mN, mNother, mSsParams, mPrng.get());
    sspmtReceiver.setTimer(getTimer());
    co_await sspmtReceiver.setup(chl);

    const u64 nSlots = sspmtReceiver.getLayoutSize();
    co_await genRotReceiver(otReceiver, mPrng, chl, nSlots, mRot);
    setTimePoint("PsiSumReceiver::Setup (transfer OTs)");

    mSetupDone = true;
};

Proto PsiSumReceiver::run(const std::vector<oc::block> &X,
                          oc::u64 &payloadSum, Socket &chl)
{
    co_await setup(chl);

    BitVector sspmt;
    co_await sspmtReceiver.run(X, sspmt, chl);

    std::vector<u32> shares;
    co_await weightedSumReceiver(mRot, chl, sspmt, 1, shares);

    std::vector<u32> senderShares;
    co_await chl.recvResize(senderShares);

    payloadSum = (u32)(shares[0] - senderShares[0]);
    setTimePoint("Receiver::Sum");
};

Proto PsiThresholdSender::setup(Socket &chl)
{
    if (mSetupDone) {
        co_return;
    }
    sspmtSender.init(mN, mNother, mSsParams, mPrng.get());
    sspmtSender.setTimer(getTimer());
    co_await sspmtSender.setup(chl);

    const u64 nSlots = sspmtSender.getLayoutSize();
    co_await genRotSender(otSender, mPrng, chl, nSlots, mRot);
    setTimePoint("PsiThresholdSender::Setup (transfer OTs)");

    auto cir = sumThresholdCircuit(kShareBits);
    mGmw.setTimer(getTimer());
    mGmw.init(kThresholdInstances, cir, 1, 1ull << 16, 0, mPrng.get());
    co_await mGmw.generateTriple(chl);
    setTimePoint("PsiThresholdSender::Setup (comparison triples)");

    mSetupDone = true;
};

Proto PsiThresholdSender::run(const std::vector<oc::block> &Y, oc::u32 threshold,
                              Socket &chl)
{

    if (threshold == 0) {
        throw RTE_LOC;
    }

    co_await setup(chl);

    BitVector sspmt;
    co_await sspmtSender.run(Y, sspmt, chl);

    auto comm = chl.bytesSent() + chl.bytesReceived();

    std::vector<std::vector<u32>> weights(1, std::vector<u32>(sspmt.size(), 1));
    std::vector<u32> shares;
    co_await weightedSumSender(mRot, chl, sspmt, weights, shares);

    const u32 negShare = (u32)(0u - shares[0]);
    oc::Matrix<u32> negShareIn(kThresholdInstances, 1);
    oc::Matrix<u32> thresholdIn(kThresholdInstances, 1);
    for (u64 i = 0; i < kThresholdInstances; ++i) {
        negShareIn(i, 0) = negShare;
        thresholdIn(i, 0) = threshold - 1;
    }

    mGmw.setZeroInput(0);
    mGmw.setInput(1, negShareIn);
    mGmw.setInput(2, thresholdIn);
    co_await mGmw.run(chl);

    auto out = mGmw.getOutputView(0);
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

    const u64 nSlots = sspmtReceiver.getLayoutSize();
    co_await genRotReceiver(otReceiver, mPrng, chl, nSlots, mRot);
    setTimePoint("PsiThresholdReceiver::Setup (transfer OTs)");

    auto cir = sumThresholdCircuit(kShareBits);
    mGmw.setTimer(getTimer());
    mGmw.init(kThresholdInstances, cir, 1, 1ull << 16, 1, mPrng.get());
    co_await mGmw.generateTriple(chl);
    setTimePoint("PsiThresholdReceiver::Setup (comparison triples)");

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
    co_await weightedSumReceiver(mRot, chl, sspmt, 1, shares);

    oc::Matrix<u32> shareIn(kThresholdInstances, 1);
    for (u64 i = 0; i < kThresholdInstances; ++i) {
        shareIn(i, 0) = shares[0];
    }

    mGmw.setInput(0, shareIn);
    mGmw.setZeroInput(1);
    mGmw.setZeroInput(2);
    co_await mGmw.run(chl);

    auto out = mGmw.getOutputView(0);
    u8 senderShare;
    co_await chl.recv(senderShare);

    aboveThreshold = ((out(0, 0) ^ senderShare) & 1) != 0;
    setTimePoint("Receiver::Threshold");
};
}
