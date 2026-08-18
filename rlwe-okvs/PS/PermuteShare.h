// Implementation of the correlated oblivious permutation of Han et al.,
// "Concretely Efficient Correlated Oblivious Permutation" (ACM AsiaCCS 2026),
// https://doi.org/10.1145/3779208.3805966.
#pragma once

#include "cryptoTools/Common/BitVector.h"
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Common/Timer.h"
#include "coproto/coproto.h"
#include "libOTe/TwoChooseOne/Silent/SilentOtExtReceiver.h"
#include "libOTe/TwoChooseOne/Silent/SilentOtExtSender.h"

#include "benes.h"

namespace rlweOkvs
{
    using Proto = coproto::task<>;
    using Socket = coproto::Socket;

    oc::u64 permuteShareSwitchCount(oc::u64 n);

    class PermuteShareSender : public oc::TimerAdapter
    {
        oc::u64 mN = 0, mLogN = 0, mNumColumns = 0, mNumSwitches = 0;
        Benes mBenes;

        std::vector<oc::u8> mRotBits;
        oc::BitVector mRotChoices;
        oc::BitVector mCorrShare;

    public:

        void init(oc::u64 n);

        Proto setup(oc::SilentOtExtReceiver& ot, oc::PRNG& prng, Socket& chl);

        void setPermutation(std::vector<int> perm);

        Proto run(oc::BitVector& share, Socket& chl);

        Proto correlate(oc::SilentOtExtReceiver& ot, oc::PRNG& prng,
            Socket& chl);

        Proto apply(oc::BitVector& share, Socket& chl);

        std::vector<int> getPerm() { return mBenes.getPerm(); };
        const std::vector<int>& getPermRef() { return mBenes.getPermRef(); };

        oc::u8 firstSwitch() { return mBenes.getSwitchesAsBitVec()[0]; };
        std::vector<int> getInvPerm() { return mBenes.getInvPerm(); };
        const std::vector<int>& getInvPermRef() { return mBenes.getInvPermRef(); };
    };

    class PermuteShareReceiver : public oc::TimerAdapter
    {
        oc::u64 mN = 0, mLogN = 0, mNumColumns = 0, mNumSwitches = 0;

        std::vector<std::array<oc::u8, 2>> mSotMsgs;
        std::vector<oc::BitVector> mScratchBottom, mScratchTop;
        oc::BitVector mCorrR;
        oc::BitVector mCorrShare;

        void prepareCorrection(oc::u64 depth, oc::u64 permIdx,
            oc::BitVector& src,
            std::vector<std::array<oc::u8, 2>>& otMsgs,
            std::vector<oc::u8>& corrections);

    public:
        void init(oc::u64 n);

        Proto setup(oc::SilentOtExtSender& ot, oc::PRNG& prng, Socket& chl);

        Proto run(const oc::BitVector& inputs, oc::BitVector& outputs,
            oc::PRNG& prng, Socket& chl);

        Proto correlate(oc::SilentOtExtSender& ot, oc::PRNG& prng, Socket& chl);

        Proto apply(const oc::BitVector& inputs, oc::BitVector& outputs,
            Socket& chl);
    };
}
