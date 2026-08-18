// Adapted from the GMW implementation in https://github.com/ladnir/volepsi.
#pragma once

#include "Defines.h"
#include <vector>
#include <libOTe/TwoChooseOne/Silent/SilentOtExtSender.h>
#include <libOTe/TwoChooseOne/Silent/SilentOtExtReceiver.h>

namespace volePSI
{
    class SilentTripleGen : public oc::TimerAdapter
    {
    public:
        bool mHasBase = false;
        u64 mN, mBatchSize, mNumPer;
        Mode mMode;
        oc::PRNG mPrng;

        std::unique_ptr<oc::SilentOtExtSender[]> mBackingSenderOT;
        std::unique_ptr<oc::SilentOtExtReceiver[]> mBackingRecverOT;

        span<oc::SilentOtExtSender> mSenderOT;
        span<oc::SilentOtExtReceiver> mRecverOT;

        void init(u64 n, u64 batchSize, u64 numThreads, Mode mode, block seed);

        Proto expand(coproto::Socket& chl);

        bool hasBaseOts()  { return mHasBase; }

        Proto generateBaseCors(oc::PRNG& prng, coproto::Socket& chl);

        span<block> mMult, mAdd;
        std::vector<block> mAVec, mBVec, mDVec;
        oc::BitVector mCBitVec;
    };
}
