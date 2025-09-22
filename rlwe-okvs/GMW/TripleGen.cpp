#include "TripleGen.h"
#include "cryptoTools/Network/IOService.h"
#include "cryptoTools/Network/Session.h"

#include "libOTe/Base/BaseOT.h"

#include "libOTe/TwoChooseOne/Iknp/IknpOtExtReceiver.h"
#include "libOTe/TwoChooseOne/Iknp/IknpOtExtSender.h"

using namespace oc;

namespace volePSI
{
	void TripleGen::init(u64 n, u64 batchSize, u64 numThreads, Mode mode, block seed)
	{
		mBatchSize = batchSize;

		mNumBatchs = (n + mBatchSize - 1) / mBatchSize;
		mNumPer = oc::roundUpTo((n + mNumBatchs - 1) / mNumBatchs, 128);

		mN = mNumBatchs * mNumPer;
		mMode = mode;
		mPrng.SetSeed(seed);
        
		auto noise = SdNoiseDistribution::Stationary;
		auto sectype = SilentSecType::SemiHonest;

		if (mode & Mode::Sender)
		{
			mSenderOT.configure(mNumPer, 2, numThreads, sectype, noise);
		}

		if (mode & Mode::Receiver)
		{
			mRecverOT.configure(mNumPer, 2, numThreads, sectype, noise);
		}			
	}

#ifndef ENABLE_SSE

	// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_shuffle_epi8&ig_expand=1038,900,6922,6328
	inline block _mm_shuffle_epi8(const block& a, const block& b)
	{
		block ret;
		for (u64 i = 0; i < 16; ++i)
		{
			u8 bb = b.get<u8>()[i];
			if (bb & 128)
				ret.set<u8>(i, 0);
			else
			{
				u8 idx = bb & 15;
				ret.set<u8>(i, idx);
			}
		}

		return ret;
	}

	// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_slli_epi16&ig_expand=1038,900,6922,6328,6470
	inline block _mm_slli_epi16(const block& a, int imm)
	{
		block ret;
		for (u64 i = 0; i < 8; ++i)
		{
			ret.set<u16>(i, a.get<u16>(i) << imm);
		}

		return ret;
	}

	// https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#text=_mm_movemask_epi8&ig_expand=1038,900,6922,6328,6470,4836
	inline int _mm_movemask_epi8(const block& a)
	{
		int ret = 0;
		for (u64 i = 0; i < 16; ++i)
		{
			ret |= int(a.get<u8>()[i] >> 7) << i;
		}

		return ret;
	}

#endif

	Proto TripleGen::run(coproto::Socket& chl)
	{
		auto aIter16 = (u16*)nullptr;
		auto bIter16 = (u16*)nullptr;
		auto aIter = oc::BitIterator{};
		auto bIter = oc::BitIterator{};
		auto j = u64{};
		auto sendMsg = std::vector<std::array<block, 2>>{};
		auto recvMsg = std::vector<block>{};
		auto recvOtChoiceBits = oc::BitVector{};
		auto shuffle = std::array<block, 16>{};

		memset(shuffle.data(), 1 << 7, sizeof(*shuffle.data()) * shuffle.size());
		for (u64 i = 0; i < 16; ++i)
			shuffle[i].set<u8>(i, 0);

		setTimePoint("TripleGen::begin");
		if (mMode & Mode::Sender)
		{
			mAVec.resize(mN / 128);
			mBVec.resize(mN / 128);
			mMult = mAVec;
			mAdd = mBVec;

			aIter16 = (u16*)mAVec.data();
			bIter16 = (u16*)mBVec.data();
			aIter = oc::BitIterator((u8*)mAVec.data(), 0);
			bIter = oc::BitIterator((u8*)mBVec.data(), 0);

			assert(mNumPer % 16 == 0);
			sendMsg.resize(mNumPer);
			for (j = 0; j < mNumBatchs; ++j)
			{
                co_await mSenderOT.genBaseCors({}, mPrng, chl);
				co_await mSenderOT.silentSend(sendMsg, mPrng, chl);

				for (u64 i = 0; i < sendMsg.size(); i += 16)
				{					
					block a00 = _mm_shuffle_epi8(sendMsg[i + 0][0], shuffle[0]);
					block a01 = _mm_shuffle_epi8(sendMsg[i + 1][0], shuffle[1]);
					block a02 = _mm_shuffle_epi8(sendMsg[i + 2][0], shuffle[2]);
					block a03 = _mm_shuffle_epi8(sendMsg[i + 3][0], shuffle[3]);
					block a04 = _mm_shuffle_epi8(sendMsg[i + 4][0], shuffle[4]);
					block a05 = _mm_shuffle_epi8(sendMsg[i + 5][0], shuffle[5]);
					block a06 = _mm_shuffle_epi8(sendMsg[i + 6][0], shuffle[6]);
					block a07 = _mm_shuffle_epi8(sendMsg[i + 7][0], shuffle[7]);
					block a08 = _mm_shuffle_epi8(sendMsg[i + 8][0], shuffle[8]);
					block a09 = _mm_shuffle_epi8(sendMsg[i + 9][0], shuffle[9]);
					block a10 = _mm_shuffle_epi8(sendMsg[i + 10][0], shuffle[10]);
					block a11 = _mm_shuffle_epi8(sendMsg[i + 11][0], shuffle[11]);
					block a12 = _mm_shuffle_epi8(sendMsg[i + 12][0], shuffle[12]);
					block a13 = _mm_shuffle_epi8(sendMsg[i + 13][0], shuffle[13]);
					block a14 = _mm_shuffle_epi8(sendMsg[i + 14][0], shuffle[14]);
					block a15 = _mm_shuffle_epi8(sendMsg[i + 15][0], shuffle[15]);

					block b00 = _mm_shuffle_epi8(sendMsg[i + 0][1], shuffle[0]);
					block b01 = _mm_shuffle_epi8(sendMsg[i + 1][1], shuffle[1]);
					block b02 = _mm_shuffle_epi8(sendMsg[i + 2][1], shuffle[2]);
					block b03 = _mm_shuffle_epi8(sendMsg[i + 3][1], shuffle[3]);
					block b04 = _mm_shuffle_epi8(sendMsg[i + 4][1], shuffle[4]);
					block b05 = _mm_shuffle_epi8(sendMsg[i + 5][1], shuffle[5]);
					block b06 = _mm_shuffle_epi8(sendMsg[i + 6][1], shuffle[6]);
					block b07 = _mm_shuffle_epi8(sendMsg[i + 7][1], shuffle[7]);
					block b08 = _mm_shuffle_epi8(sendMsg[i + 8][1], shuffle[8]);
					block b09 = _mm_shuffle_epi8(sendMsg[i + 9][1], shuffle[9]);
					block b10 = _mm_shuffle_epi8(sendMsg[i + 10][1], shuffle[10]);
					block b11 = _mm_shuffle_epi8(sendMsg[i + 11][1], shuffle[11]);
					block b12 = _mm_shuffle_epi8(sendMsg[i + 12][1], shuffle[12]);
					block b13 = _mm_shuffle_epi8(sendMsg[i + 13][1], shuffle[13]);
					block b14 = _mm_shuffle_epi8(sendMsg[i + 14][1], shuffle[14]);
					block b15 = _mm_shuffle_epi8(sendMsg[i + 15][1], shuffle[15]);

					a00 = a00 ^ a08;
					a01 = a01 ^ a09;
					a02 = a02 ^ a10;
					a03 = a03 ^ a11;
					a04 = a04 ^ a12;
					a05 = a05 ^ a13;
					a06 = a06 ^ a14;
					a07 = a07 ^ a15;

					b00 = b00 ^ b08;
					b01 = b01 ^ b09;
					b02 = b02 ^ b10;
					b03 = b03 ^ b11;
					b04 = b04 ^ b12;
					b05 = b05 ^ b13;
					b06 = b06 ^ b14;
					b07 = b07 ^ b15;

					a00 = a00 ^ a04;
					a01 = a01 ^ a05;
					a02 = a02 ^ a06;
					a03 = a03 ^ a07;

					b00 = b00 ^ b04;
					b01 = b01 ^ b05;
					b02 = b02 ^ b06;
					b03 = b03 ^ b07;

					a00 = a00 ^ a02;
					a01 = a01 ^ a03;

					b00 = b00 ^ b02;
					b01 = b01 ^ b03;

					a00 = a00 ^ a01;
					b00 = b00 ^ b01;

					a00 = _mm_slli_epi16(a00, 7);
					b00 = _mm_slli_epi16(b00, 7);

					u16 ap = _mm_movemask_epi8(a00);
					u16 bp = _mm_movemask_epi8(b00);

					*aIter16++ = ap ^ bp;
					*bIter16++ = ap;
				}
			}
			sendMsg = {};
		}
		else
		{
			mDVec.resize(mN / 128);
			mCBitVec.resize(0);
			mCBitVec.reserve(mN);
			mAdd = mDVec;
			mMult = span<block>((block*)mCBitVec.data(), mAdd.size());
			recvMsg.resize(mNumPer);
			recvOtChoiceBits.resize(mNumPer);
			aIter = oc::BitIterator((u8*)mAdd.data(), 0);
			aIter16 = (u16*)mAdd.data();

			for (j = 0; j < mNumBatchs; ++j)
			{
                co_await mRecverOT.genBaseCors(mPrng, chl);
				co_await mRecverOT.silentReceive(recvOtChoiceBits, recvMsg, mPrng, chl);

				mCBitVec.append(recvOtChoiceBits);


				for (u64 i = 0; i < recvMsg.size(); i += 16)
				{
					block a00 = _mm_shuffle_epi8(recvMsg[i + 0], shuffle[0]);
					block a01 = _mm_shuffle_epi8(recvMsg[i + 1], shuffle[1]);
					block a02 = _mm_shuffle_epi8(recvMsg[i + 2], shuffle[2]);
					block a03 = _mm_shuffle_epi8(recvMsg[i + 3], shuffle[3]);
					block a04 = _mm_shuffle_epi8(recvMsg[i + 4], shuffle[4]);
					block a05 = _mm_shuffle_epi8(recvMsg[i + 5], shuffle[5]);
					block a06 = _mm_shuffle_epi8(recvMsg[i + 6], shuffle[6]);
					block a07 = _mm_shuffle_epi8(recvMsg[i + 7], shuffle[7]);
					block a08 = _mm_shuffle_epi8(recvMsg[i + 8], shuffle[8]);
					block a09 = _mm_shuffle_epi8(recvMsg[i + 9], shuffle[9]);
					block a10 = _mm_shuffle_epi8(recvMsg[i + 10], shuffle[10]);
					block a11 = _mm_shuffle_epi8(recvMsg[i + 11], shuffle[11]);
					block a12 = _mm_shuffle_epi8(recvMsg[i + 12], shuffle[12]);
					block a13 = _mm_shuffle_epi8(recvMsg[i + 13], shuffle[13]);
					block a14 = _mm_shuffle_epi8(recvMsg[i + 14], shuffle[14]);
					block a15 = _mm_shuffle_epi8(recvMsg[i + 15], shuffle[15]);

					a00 = a00 ^ a08;
					a01 = a01 ^ a09;
					a02 = a02 ^ a10;
					a03 = a03 ^ a11;
					a04 = a04 ^ a12;
					a05 = a05 ^ a13;
					a06 = a06 ^ a14;
					a07 = a07 ^ a15;

					a00 = a00 ^ a04;
					a01 = a01 ^ a05;
					a02 = a02 ^ a06;
					a03 = a03 ^ a07;

					a00 = a00 ^ a02;
					a01 = a01 ^ a03;

					a00 = a00 ^ a01;

					a00 = _mm_slli_epi16(a00, 7);

					u16 ap = _mm_movemask_epi8(a00);

					*aIter16++ = ap;
				}
			}
		}
		setTimePoint("TripleGen::end");
	}
}
