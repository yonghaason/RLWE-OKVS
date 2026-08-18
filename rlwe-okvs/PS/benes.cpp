#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>
#include <stack>
#include <fstream>
#include <numeric>
#include <random>
#include <cryptoTools/Common/Defines.h>
#include <cryptoTools/Common/BitVector.h>
#include "benes.h"

using namespace std;
using namespace oc;

namespace rlweOkvs {

void Benes::initialize(int N, vector<int> perm)
{
	mN = N;
	mLogN = int(ceil(log2(N)));
	mNumColumns = 2 * mLogN - 1;

	if (perm.size())
	{
		assert(perm.size() == N);
		mPerm = perm;
	}
	else
	{
		mPerm.resize(N);
		std::iota(mPerm.begin(), mPerm.end(), 0);

		std::random_device rd;
		std::shuffle(mPerm.begin(), mPerm.end(), std::default_random_engine(rd()));
	}

	mInvPerm.resize(N);
	for (size_t i = 0; i < N; i++)
	{
		mInvPerm[mPerm[i]] = i;
	}

	permInner.resize(N);
	invPermInner.resize(N);

	mSwitches.resize(mNumColumns);
	for (int i = 0; i < mNumColumns; ++i)
		mSwitches[i].resize(mN / 2);
}

void Benes::DFS(int idx, int route,
								vector<char> &path)
{
	stack<pair<int, int>> st;
	st.push({idx, route});
	pair<int, int> pr;
	while (!st.empty())
	{
		pr = st.top();
		st.pop();
		path[pr.first] = pr.second;
		if (path[pr.first ^ 1] < 0)
			st.push({pr.first ^ 1,
							 pr.second ^ 1});

		idx = permInner[invPermInner[pr.first] ^ 1];
		if (path[idx] < 0)
			st.push({idx, pr.second ^ 1});
	}
}

int shuffle(int i, int n)
{
	return ((i & 1) << (n - 1)) | (i >> 1);
}

void Benes::genBenesRouteInner(int depth, int permIdx, const vector<int> &src,
															 const vector<int> &dest)
{
	int i, j, x, s;
	vector<int> &bottom1 = mRouteBottom1[depth];
	vector<int> &top1 = mRouteTop1[depth];
	bottom1.clear();
	top1.clear();
	int subNetSize = src.size();

	int coDepth = mLogN - depth;

	if (subNetSize == 2)
	{
		if (coDepth == 1)
			mSwitches[depth][permIdx] = src[0] != dest[0];
		else
			mSwitches[depth + 1][permIdx] = src[0] != dest[0];
		return;
	}

	if (subNetSize == 3)
	{
		if (src[0] == dest[0])
		{
			mSwitches[depth][permIdx] = 0;
			mSwitches[depth + 2][permIdx] = 0;
			if (src[1] == dest[1])
				mSwitches[depth + 1][permIdx] = 0;
			else
				mSwitches[depth + 1][permIdx] = 1;
		}

		if (src[0] == dest[1])
		{
			mSwitches[depth][permIdx] = 0;
			mSwitches[depth + 2][permIdx] = 1;
			if (src[1] == dest[0])
				mSwitches[depth + 1][permIdx] = 0;
			else
				mSwitches[depth + 1][permIdx] = 1;
		}

		if (src[0] == dest[2])
		{
			mSwitches[depth][permIdx] = 1;
			mSwitches[depth + 1][permIdx] = 1;
			if (src[1] == dest[0])
				mSwitches[depth + 2][permIdx] = 0;
			else
				mSwitches[depth + 2][permIdx] = 1;
		}
		return;
	}

	int levels = 2 * coDepth - 1;

	vector<int> &bottom2 = mRouteBottom2[depth];
	vector<int> &top2 = mRouteTop2[depth];
	bottom2.assign(subNetSize / 2, 0);
	top2.assign(int(ceil(subNetSize * 0.5)), 0);

	for (i = 0; i < subNetSize; ++i)
		invPermInner[src[i]] = i;

	for (i = 0; i < subNetSize; ++i)
		permInner[i] = invPermInner[dest[i]];

	for (i = 0; i < subNetSize; ++i)
		invPermInner[permInner[i]] = i;

	vector<char> path(subNetSize, -1);
	if (subNetSize % 2 == 1)
	{
		path[subNetSize - 1] = 1;
		path[permInner[subNetSize - 1]] = 1;
		if (permInner[subNetSize - 1] != subNetSize - 1)
		{
			int idx = permInner[invPermInner[subNetSize - 1] ^ 1];
			DFS(idx, 0, path);
		}
	}

	for (i = 0; i < subNetSize; ++i)
	{
		if (path[i] < 0)
		{
			DFS(i, 0, path);
		}
	}

	bottom1.resize(subNetSize / 2);
	top1.resize(subNetSize / 2 + (subNetSize % 2));
	for (i = 0; i < subNetSize - 1; i += 2)
	{
		s = path[i];
		mSwitches[depth][permIdx + i / 2] = s;
		bottom1[i / 2] = src[i + s];
		top1[i / 2] = src[i + 1 - s];
	}
	if (subNetSize % 2 == 1)
	{
		top1[subNetSize / 2] = src[subNetSize - 1];
	}

	for (i = 0; i < subNetSize - 1; i += 2)
	{
		s = mSwitches[depth + levels - 1][permIdx + i / 2] = path[permInner[i]];
		bottom2[i / 2] = src[permInner[i + s]];
		top2[i / 2] = src[permInner[i + 1 - s]];
	}

	int idx = int(ceil(subNetSize * 0.5));
	if (subNetSize % 2 == 1)
	{
		top2[idx - 1] = dest[subNetSize - 1];
	}

	genBenesRouteInner(depth + 1, permIdx, bottom1, bottom2);
	genBenesRouteInner(depth + 1, permIdx + subNetSize / 4, top1, top2);
}

void Benes::benesEval(vector<int> &src, int depth, int permIdx)
{
	int i, j, x, s;
	vector<int> bottom1;
	vector<int> top1;
	int subNetSize = src.size();
	int coDepth = mLogN - depth;
	int temp;

	if (subNetSize == 2)
	{
		if (coDepth == 1)
		{
			if (mSwitches[depth][permIdx] == 1)
			{
				temp = src[0];
				src[0] = src[1];
				src[1] = temp;
			}
		}
		else if (mSwitches[depth + 1][permIdx] == 1)
		{
			temp = src[0];
			src[0] = src[1];
			src[1] = temp;
		}
		return;
	}

	if (subNetSize == 3)
	{
		if (mSwitches[depth][permIdx] == 1)
		{
			temp = src[0];
			src[0] = src[1];
			src[1] = temp;
		}
		if (mSwitches[depth + 1][permIdx] == 1)
		{
			temp = src[1];
			src[1] = src[2];
			src[2] = temp;
		}
		if (mSwitches[depth + 2][permIdx] == 1)
		{
			temp = src[0];
			src[0] = src[1];
			src[1] = temp;
		}
		return;
	}

	int levels = 2 * coDepth - 1;

	for (i = 0; i < subNetSize - 1; i += 2)
	{
		int s = mSwitches[depth][permIdx + i / 2];
		for (j = 0; j < 2; ++j)
		{
			x = shuffle((i | j) ^ s, coDepth);
			if (x < subNetSize / 2)
				bottom1.push_back(src[i | j]);
			else
				top1.push_back(src[i | j]);
		}
	}
	if (subNetSize % 2 == 1)
	{
		top1.push_back(src[subNetSize - 1]);
	}

	benesEval(bottom1, depth + 1, permIdx);
	benesEval(top1, depth + 1, permIdx + subNetSize / 4);

	for (i = 0; i < subNetSize - 1; i += 2)
	{
		s = mSwitches[depth + levels - 1][permIdx + i / 2];
		for (j = 0; j < 2; ++j)
		{
			x = shuffle((i | j) ^ s, coDepth);
			if (x < subNetSize / 2)
				src[i | j] = bottom1[x];
			else
			{
				src[i | j] = top1[i / 2];
			}
		}
	}

	int idx = int(ceil(subNetSize * 0.5));
	if (subNetSize % 2 == 1)
	{
		src[subNetSize - 1] = top1[idx - 1];
	}
}

void Benes::reserveScratch()
{
	mRouteBottom1.resize(mLogN + 2);
	mRouteTop1.resize(mLogN + 2);
	mRouteBottom2.resize(mLogN + 2);
	mRouteTop2.resize(mLogN + 2);
	mEvalBottom.resize(mLogN + 2);
	mEvalTop.resize(mLogN + 2);
	for (int i = 0; i <= mLogN + 1; ++i)
	{
		mRouteBottom1[i].reserve(mN);
		mRouteTop1[i].reserve(mN);
		mRouteBottom2[i].reserve(mN);
		mRouteTop2[i].reserve(mN);
		mEvalBottom[i].reserve(mN);
		mEvalTop[i].reserve(mN);
	}
}

void Benes::genBenesRoute()
{
	reserveScratch();
	vector<int> src(mN);
	std::iota(src.begin(), src.end(), 0);
	genBenesRouteInner(0, 0, src, mPerm);
}

void Benes::benesMaskedEval(oc::BitVector &src,
														std::vector<oc::u8> &otMsgs,
														int depth, int permIdx)
{
	int levels, i, j, x, s;
	oc::BitVector &bottom1 = mEvalBottom[depth];
	oc::BitVector &top1 = mEvalTop[depth];
	bottom1.resize(0);
	top1.resize(0);
	int subNetSize = src.size();
	int coDepth = mLogN - depth;
	oc::u8 temp, temp_int;

	if (subNetSize == 2)
	{
		if (coDepth == 1)
		{
			temp_int = otMsgs[(depth) * (mN / 2) + (permIdx)];
			src[0] = src[0] ^ temp_int;
			src[1] = src[1] ^ temp_int;
			if (mSwitches[depth][permIdx] == 1)
			{
				temp = src[0];
				src[0] = src[1];
				src[1] = temp;
			}
		}
		else
		{
			temp_int = otMsgs[(depth + 1) * (mN / 2) + (permIdx)];
			src[0] = src[0] ^ temp_int;
			src[1] = src[1] ^ temp_int;
			if (mSwitches[depth + 1][permIdx] == 1)
			{
				temp = src[0];
				src[0] = src[1];
				src[1] = temp;
			}
		}
		return;
	}

	if (subNetSize == 3)
	{
		temp_int = otMsgs[(depth) * (mN / 2) + (permIdx)];
		src[0] = src[0] ^ temp_int;
		src[1] = src[1] ^ temp_int;
		if (mSwitches[depth][permIdx] == 1)
		{
			temp = src[0];
			src[0] = src[1];
			src[1] = temp;
		}

		temp_int = otMsgs[(depth + 1) * (mN / 2) + (permIdx)];
		src[1] = src[1] ^ temp_int;
		src[2] = src[2] ^ temp_int;
		if (mSwitches[depth + 1][permIdx] == 1)
		{
			temp = src[1];
			src[1] = src[2];
			src[2] = temp;
		}

		temp_int = otMsgs[(depth + 2) * (mN / 2) + (permIdx)];
		src[0] = src[0] ^ temp_int;
		src[1] = src[1] ^ temp_int;
		if (mSwitches[depth + 2][permIdx] == 1)
		{
			temp = src[0];
			src[0] = src[1];
			src[1] = temp;
		}

		return;
	}

	levels = 2 * coDepth - 1;

	bottom1.resize(subNetSize / 2);
	top1.resize(subNetSize / 2 + (subNetSize % 2));
	for (i = 0; i < subNetSize - 1; i += 2)
	{
		s = mSwitches[depth][permIdx + i / 2];
		temp_int = otMsgs[(depth) * (mN / 2) + (permIdx + i / 2)];
		src[i] = src[i] ^ temp_int;
		src[i ^ 1] = src[i ^ 1] ^ temp_int;
		bottom1[i / 2] = (oc::u8)src[i + s];
		top1[i / 2] = (oc::u8)src[i + 1 - s];
	}
	if (subNetSize % 2 == 1)
	{
		top1[subNetSize / 2] = (oc::u8)src[subNetSize - 1];
	}

	benesMaskedEval(bottom1, otMsgs, depth + 1, permIdx);
	benesMaskedEval(top1, otMsgs, depth + 1, permIdx + subNetSize / 4);

	for (i = 0; i < subNetSize - 1; i += 2)
	{
		s = mSwitches[depth + levels - 1][permIdx + i / 2];

		src[i + s] = (oc::u8)bottom1[i / 2];
		src[i + 1 - s] = (oc::u8)top1[i / 2];

		temp_int = otMsgs[(depth + levels - 1) * (mN / 2) + (permIdx + i / 2)];
		src[i] = src[i] ^ temp_int;
		src[i ^ 1] = src[i ^ 1] ^ temp_int;
	}

	int idx = int(ceil(subNetSize * 0.5));
	if (subNetSize % 2 == 1)
	{
		src[subNetSize - 1] = top1[idx - 1];
	}
}

oc::BitVector Benes::getSwitchesAsBitVec()
{
	oc::BitVector stretchedSwitches(mNumColumns * (mN / 2));
	oc::u8 *out = stretchedSwitches.data();
	for (int j = 0; j < mNumColumns; ++j)
	{
		const size_t base = (size_t)(mN / 2) * j;
		for (int i = 0; i < mN / 2; ++i)
		{
			const size_t k = base + i;
			out[k >> 3] |= (oc::u8)((mSwitches[j][i] & 1) << (k & 7));
		}
	}
	return stretchedSwitches;
}
}
