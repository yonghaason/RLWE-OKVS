#pragma once
#include "cryptoTools/Common/Defines.h"
#include "cryptoTools/Common/block.h"
#include "cryptoTools/Common/BitVector.h"
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Common/Timer.h"
#include "cryptoTools/Network/Channel.h"
#include "coproto/coproto.h"
#include "rlwe-okvs/sspmt.h"
#include "PS/PermuteShare.h"
#include "libOTe/TwoChooseOne/Silent/SilentOtExtReceiver.h"
#include "libOTe/TwoChooseOne/Silent/SilentOtExtSender.h"
#include "seal/seal.h"

namespace rlweOkvs
{
    using Proto = coproto::task<>;
    using Socket = coproto::Socket;

    struct RotStore
    {
        std::vector<std::array<oc::block, 2>> mSend;
        std::vector<oc::block> mRecv;
        oc::BitVector mChoices;
    };

    Proto genRotSender(oc::SilentOtExtSender& ot, oc::PRNG& prng, Socket& chl,
        oc::u64 count, RotStore& store);
    Proto genRotReceiver(oc::SilentOtExtReceiver& ot, oc::PRNG& prng,
        Socket& chl, oc::u64 count, RotStore& store);

    Proto weightedSumSender(
        const RotStore& rot, Socket& chl,
        const oc::BitVector& shareBits,
        const std::vector<std::vector<oc::u32>>& weights,
        std::vector<oc::u32>& shares);

    Proto weightedSumReceiver(
        const RotStore& rot, Socket& chl,
        const oc::BitVector& shareBits, oc::u64 numWeights,
        std::vector<oc::u32>& shares);

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

    class PsuSender : public PsoBase
    {
        SspmtSender sspmtSender;
        oc::SilentOtExtSender otSender;
        oc::SilentOtExtReceiver otSwitchReceiver;
        PermuteShareSender psSender;

    public:

        Proto setup(Socket& chl);

        Proto run(const std::vector<oc::block>& Y, Socket& chl);
    };

    class PsuReceiver : public PsoBase
    {
        SspmtReceiver sspmtReceiver;
        oc::SilentOtExtReceiver otReceiver;
        oc::SilentOtExtSender otSwitchSender;
        PermuteShareReceiver psReceiver;

    public:

        Proto setup(Socket& chl);

        Proto run(const std::vector<oc::block>& X,
            std::vector<oc::block>& D, Socket& chl);
    };

    class PsiCardSender : public PsoBase
    {
        SspmtSender sspmtSender;
        oc::SilentOtExtSender otSender;

    public:

        Proto setup(Socket& chl);

        Proto run(const std::vector<oc::block>& Y, Socket& chl);
    };

    class PsiCardReceiver : public PsoBase
    {
        SspmtReceiver sspmtReceiver;
        oc::SilentOtExtReceiver otReceiver;

    public:

        Proto setup(Socket& chl);

        Proto run(const std::vector<oc::block>& X,
            oc::u64& cardinality, Socket& chl);
    };

    class PsiSumSender : public PsoBase
    {
        SspmtSender sspmtSender;
        oc::SilentOtExtSender otSender;

    public:

        Proto setup(Socket& chl);

        Proto run(const std::vector<oc::block>& Y,
            const std::vector<oc::u32>& payloads, Socket& chl);
    };

    class PsiSumReceiver : public PsoBase
    {
        SspmtReceiver sspmtReceiver;
        oc::SilentOtExtReceiver otReceiver;

    public:

        Proto setup(Socket& chl);

        Proto run(const std::vector<oc::block>& X,
            oc::u64& payloadSum, Socket& chl);
    };

    class PsiThresholdSender : public PsoBase
    {
        SspmtSender sspmtSender;
        oc::SilentOtExtSender otSender;

        volePSI::Gmw mGmw;

    public:

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

        Proto setup(Socket& chl);

        Proto run(const std::vector<oc::block>& X, oc::u32 threshold,
            bool& aboveThreshold, Socket& chl);
    };
};
