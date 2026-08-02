#pragma once

#include <vector>

#include <cryptoTools/Common/BitVector.h>

namespace rlweOkvs {

class Benes
{
	int mN;
	int mNumColumns;
	int mLogN;

	std::vector<std::vector<char>> mSwitches;

	std::vector<int> mPerm;
	std::vector<int> mInvPerm;

	std::vector<int> permInner;
	std::vector<int> invPermInner;

	void DFS(int idx, int route, std::vector<char> &path);

	void genBenesRouteInner(int depth, int permIdx, const std::vector<int> &src,
													const std::vector<int> &dest);

public:

	void initialize(int N, std::vector<int> perm = std::vector<int>());

	void genBenesRoute();

	void benesEval(std::vector<int> &vec, int depth = 0, int permIdx = 0);

	void benesMaskedEval(oc::BitVector &src,
											std::vector<std::vector<std::array<oc::u8, 2>>> &otMsgs,
											int depth = 0, int permIdx = 0);


	oc::BitVector getSwitchesAsBitVec();

	std::vector<int> getPerm() {return mPerm;};
	std::vector<int> getInvPerm() {return mInvPerm;};
};

}  // namespace rlweOkvs
