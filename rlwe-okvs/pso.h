#pragma once
#include "cryptoTools/Common/Defines.h"
#include "cryptoTools/Common/block.h"
#include "cryptoTools/Common/BitVector.h"
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Common/Timer.h"
#include "cryptoTools/Network/Channel.h"
#include "coproto/coproto.h"
#include "rlwe-okvs/sspmt.h"
#include "libOTe/TwoChooseOne/Silent/SilentOtExtReceiver.h"
#include "libOTe/TwoChooseOne/Silent/SilentOtExtSender.h"
#include "seal/seal.h"

namespace rlweOkvs
{
    using Proto = coproto::task<>;
    using Socket = coproto::Socket;

    // Random OTs held by a party, generated in the offline phase. Their count
    // is a public function of the parameters, and their choice bits are
    // random, so nothing here depends on either input; the online phase just
    // derandomizes them (one bit per OT).
    struct RotStore
    {
        std::vector<std::array<oc::block, 2>> mSend;   // sender side
        std::vector<oc::block> mRecv;                  // receiver side
        oc::BitVector mChoices;                        // receiver side
    };

    Proto genRotSender(oc::SilentOtExtSender& ot, oc::PRNG& prng, Socket& chl,
        oc::u64 count, RotStore& store);
    Proto genRotReceiver(oc::SilentOtExtReceiver& ot, oc::PRNG& prng,
        Socket& chl, oc::u64 count, RotStore& store);

    // Turn the ssPMT output -- XOR shares b_s = u_s ^ v_s of the per-slot
    // membership bits -- into additive shares of sum_s b_s * w_s (mod 2^32)
    // for each of the given weight vectors, using one random OT per (weight,
    // slot). The receiver's share minus the sender's is the weighted sum, so
    // nothing about the individual bits is revealed. With w = 1 this is the
    // intersection cardinality; with w = the payload at that slot it is the
    // payload sum over the intersection.
    //
    // The OTs come from the offline phase, so all that happens here is the
    // choice-bit correction and one u32 per OT.
    Proto weightedSumSender(
        const RotStore& rot, Socket& chl,
        const oc::BitVector& shareBits,
        const std::vector<std::vector<oc::u32>>& weights,
        std::vector<oc::u32>& shares);

    Proto weightedSumReceiver(
        const RotStore& rot, Socket& chl,
        const oc::BitVector& shareBits, oc::u64 numWeights,
        std::vector<oc::u32>& shares);

    // Common state of the ssPMT-based private set operations.
    struct PsoBase : public oc::TimerAdapter
    {
        uint32_t mN = 0;
        uint32_t mNother = 0;
        oc::PRNG mPrng;
        sspmtParams mSsParams;

        void initWithParam(uint32_t n, uint32_t nOther,
            sspmtParams ssParams, oc::block seed) {
            mN = n;
            mNother = nOther;
            mSsParams = ssParams;
            mPrng.SetSeed(seed);
        };

        void init(uint32_t n, uint32_t nOther, oc::block seed) {
            mN = n;
            mNother = nOther;
            mSsParams.initialize(n);
            mPrng.SetSeed(seed);
        };

        bool mSetupDone = false;
        RotStore mRot;
    };

    // Private set union. The sender transfers, per layout slot, either the
    // item sitting there or garbage, so the receiver ends up with Y \ X.
    //
    // NOTE: unlike the other operations below, this one exposes the layout
    // position of every transferred element. The receiver can hash its own
    // items to slot columns, so those positions -- not the membership bits --
    // are what a Chandran-style analysis attacks. Closing it needs the
    // transfer to be indexed by sender-item order rather than by slot, i.e. a
    // permute+share step on the membership shares before the OT.
    class PsuSender : public PsoBase
    {
        SspmtSender sspmtSender;
        oc::SilentOtExtSender otSender;

    public:
        // Offline phase: the ss-PMT's input-independent correlated randomness.
        // run() does it inline when skipped.
        Proto setup(Socket& chl);

        Proto run(const std::vector<oc::block>& Y, Socket& chl);
    };

    class PsuReceiver : public PsoBase
    {
        SspmtReceiver sspmtReceiver;
        oc::SilentOtExtReceiver otReceiver;

    public:
        // Offline phase: the ss-PMT's input-independent correlated randomness.
        // run() does it inline when skipped.
        Proto setup(Socket& chl);

        Proto run(const std::vector<oc::block>& X,
            std::vector<oc::block>& D, Socket& chl);
    };

    // Intersection cardinality: |X n Y| to the receiver only.
    class PsiCardSender : public PsoBase
    {
        SspmtSender sspmtSender;
        oc::SilentOtExtSender otSender;

    public:
        // Offline phase: the ss-PMT's input-independent correlated randomness.
        // run() does it inline when skipped.
        Proto setup(Socket& chl);

        Proto run(const std::vector<oc::block>& Y, Socket& chl);
    };

    class PsiCardReceiver : public PsoBase
    {
        SspmtReceiver sspmtReceiver;
        oc::SilentOtExtReceiver otReceiver;

    public:
        // Offline phase: the ss-PMT's input-independent correlated randomness.
        // run() does it inline when skipped.
        Proto setup(Socket& chl);

        Proto run(const std::vector<oc::block>& X,
            oc::u64& cardinality, Socket& chl);
    };

    // Intersection cardinality together with the sum of the sender's payloads
    // over the intersection.
    class PsiCardSumSender : public PsoBase
    {
        SspmtSender sspmtSender;
        oc::SilentOtExtSender otSender;

    public:
        // Offline phase: the ss-PMT's input-independent correlated randomness.
        // run() does it inline when skipped.
        Proto setup(Socket& chl);

        Proto run(const std::vector<oc::block>& Y,
            const std::vector<oc::u32>& payloads, Socket& chl);
    };

    class PsiCardSumReceiver : public PsoBase
    {
        SspmtReceiver sspmtReceiver;
        oc::SilentOtExtReceiver otReceiver;

    public:
        // Offline phase: the ss-PMT's input-independent correlated randomness.
        // run() does it inline when skipped.
        Proto setup(Socket& chl);

        Proto run(const std::vector<oc::block>& X,
            oc::u64& cardinality, oc::u64& payloadSum, Socket& chl);
    };

    // Threshold: only the bit 1(|X n Y| >= t), the cardinality itself stays
    // secret-shared and is compared inside a GMW circuit.
    class PsiThresholdSender : public PsoBase
    {
        SspmtSender sspmtSender;
        oc::SilentOtExtSender otSender;
        // The comparison circuit: one instance of real input, sized and
        // supplied with triples in the offline phase like everything else.
        volePSI::Gmw mGmw;

    public:
        // Offline phase: the ss-PMT's input-independent correlated randomness.
        // run() does it inline when skipped.
        Proto setup(Socket& chl);

        Proto run(const std::vector<oc::block>& Y, oc::u32 threshold,
            Socket& chl);
    };

    class PsiThresholdReceiver : public PsoBase
    {
        SspmtReceiver sspmtReceiver;
        oc::SilentOtExtReceiver otReceiver;
        volePSI::Gmw mGmw;

    public:
        // Offline phase: the ss-PMT's input-independent correlated randomness.
        // run() does it inline when skipped.
        Proto setup(Socket& chl);

        Proto run(const std::vector<oc::block>& X, oc::u32 threshold,
            bool& aboveThreshold, Socket& chl);
    };
};
