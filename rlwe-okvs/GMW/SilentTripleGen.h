#pragma once
// © 2022 Visa.
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

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

        // A * C = B + D
        span<block> mMult, mAdd;
        std::vector<block> mAVec, mBVec, mDVec;
        oc::BitVector mCBitVec;
    };
}