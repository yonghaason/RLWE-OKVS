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
// Bit width of the arithmetic shares produced by weightedSum*. Cardinalities
// and payload sums are reduced mod 2^32.
constexpr u64 kShareBits = 32;

// Number of GMW instances used for the (single) threshold comparison. The
// circuit is evaluated on one row of real input; the rest are padding so the
// silent OT backend stays in a comfortable size regime.
constexpr u64 kThresholdInstances = 128;
}  // namespace

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

    // The offline OTs have random choice bits; the receiver says which of them
    // to flip so that the effective choice is its share bit.
    BitVector flip(nSlots * k);
    co_await chl.recv(flip);

    // With the pair oriented to the receiver's choice, sending
    // delta = r0 - r1 + (1-2v)w lets it reach r0 + u(1-2v)w, whose distance to
    // the sender's r0 - vw is exactly (u ^ v) * w.
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

    // The whole permute + share is precomputed against a random permutation:
    // switch OTs, routing and both network evaluations. Online it costs one
    // blinded vector and a local permutation.
    const u64 nSlots = sspmtSender.getLayoutSize();
    psSender.setTimer(getTimer());
    psSender.init(nSlots);
    co_await psSender.correlate(otSwitchReceiver, mPrng, chl);

    // The transfer runs over item indices now, not slots, so n_y OTs suffice.
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

    // The permutation was fixed offline and at random, so instead of routing
    // one that puts our items first, we look up where it sent their slots.
    // Those positions are a uniformly random injection from the receiver's
    // point of view, so naming them reveals nothing.
    //
    // The item order still has to be randomized: the receiver does see which
    // OT index carried a new element, and that index must not be our input
    // order, or the pattern of matched items leaks by rank.
    std::vector<int> itemToSlot(mN, -1);
    for (u64 s = 0; s < nSlots; ++s) {
        if (slotToItem[s] != UINT32_MAX) itemToSlot[slotToItem[s]] = (int)s;
    }
    std::vector<u32> order(mN);
    std::iota(order.begin(), order.end(), 0u);
    for (u64 i = mN - 1; i > 0; --i) std::swap(order[i], order[mPrng.get<u64>() % (i + 1)]);

    BitVector psShare;
    co_await psSender.apply(psShare, chl);

    const auto &invPerm = psSender.getInvPermRef();
    std::vector<u32> positions(mN);
    BitVector bits(mN);
    for (u64 j = 0; j < mN; ++j) {
        const u32 slot = (u32)itemToSlot[order[j]];
        const u32 p = (u32)invPerm[slot];
        positions[j] = p;
        // psShare[p] ^ (receiver's share)[p] is the membership bit at `slot`;
        // fold in our own ss-PMT share of that slot.
        bits[j] = psShare[p] ^ sspmt[slot];
    }
    co_await chl.send(coproto::copy(positions));

    BitVector flip(mN);
    co_await chl.recv(flip);

    // Send only the message the receiver recovers when the shares agree; if
    // they disagree it xors with the other pad and gets garbage.
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

    // Permute + share, precomputed offline against a permutation only the
    // sender knows; here it is one blinded vector.
    BitVector psShare;
    co_await psReceiver.apply(sspmt, psShare, chl);

    std::vector<u32> positions;
    co_await chl.recvResize(positions);

    BitVector bits(mNother);
    for (u64 j = 0; j < mNother; ++j) bits[j] = psShare[positions[j]];

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

    // The transfer OTs are random and their count is public, so they belong
    // here too; the online phase only flips them onto the real choice bits.
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

    // The transfer OTs are random and their count is public, so they belong
    // here too; the online phase only flips them onto the real choice bits.
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

Proto PsiCardSumSender::setup(Socket &chl)
{
    if (mSetupDone) {
        co_return;
    }
    sspmtSender.init(mN, mNother, mSsParams, mPrng.get());
    sspmtSender.setTimer(getTimer());
    co_await sspmtSender.setup(chl);

    // The transfer OTs are random and their count is public, so they belong
    // here too; the online phase only flips them onto the real choice bits.
    const u64 nSlots = sspmtSender.getLayoutSize();
    co_await genRotSender(otSender, mPrng, chl, 2 * nSlots, mRot);
    setTimePoint("PsiCardSumSender::Setup (transfer OTs)");

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
    co_await weightedSumSender(mRot, chl, sspmt, weights, shares);
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

    // The transfer OTs are random and their count is public, so they belong
    // here too; the online phase only flips them onto the real choice bits.
    const u64 nSlots = sspmtReceiver.getLayoutSize();
    co_await genRotReceiver(otReceiver, mPrng, chl, 2 * nSlots, mRot);
    setTimePoint("PsiCardSumReceiver::Setup (transfer OTs)");

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
    co_await weightedSumReceiver(mRot, chl, sspmt, 2, shares);

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

    // The transfer OTs are random and their count is public, so they belong
    // here too; the online phase only flips them onto the real choice bits.
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
    co_await weightedSumSender(mRot, chl, sspmt, weights, shares);

    // The cardinality stays shared: the comparison against the threshold runs
    // inside GMW over (receiverShare - senderShare), so neither party sees the
    // count. Bundle 0 takes the receiver's share, bundle 1 the negated
    // sender's share and bundle 2 the public threshold; each party feeds its
    // own bundles and zeroes the others. The circuit tests the strict
    // "> threshold - 1", which is "at least threshold".
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
    co_await mGmw.run(chl);  // triples already generated -> online only

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

    // The transfer OTs are random and their count is public, so they belong
    // here too; the online phase only flips them onto the real choice bits.
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
    co_await mGmw.run(chl);  // triples already generated -> online only

    auto out = mGmw.getOutputView(0);
    u8 senderShare;
    co_await chl.recv(senderShare);

    aboveThreshold = ((out(0, 0) ^ senderShare) & 1) != 0;
    setTimePoint("Receiver::Threshold");
};
}  // namespace rlweOkvs
