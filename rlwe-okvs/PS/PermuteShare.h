#pragma once
// Permute + Share over single bits.
//
// F_PS: the sender holds a permutation pi over [n], the receiver holds an
// input vector x of n bits. The receiver ends with a, the sender with b, such
// that a ^ b = pi(x) -- neither learns pi(x), and the receiver learns nothing
// about pi. This is what lets the ss-PMT membership shares be handed back in
// the sender's item order without disclosing the layout position of any item,
// which is what the PSU transfer would otherwise expose.
//
// Ported from github.com/yonghaason/volepsi (volePSI/PS), restricted to bits
// and restructured so the oblivious transfers -- one per Benes switch, random
// and input-independent -- are generated in an offline phase.

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

    // Number of Benes switches for n wires: 2*ceil(log n) - 1 columns of
    // n/2 switches each. One OT per switch.
    oc::u64 permuteShareSwitchCount(oc::u64 n);

    // Holds the permutation. It is the OT *receiver*: the choice bit of switch
    // i is that switch's setting, so the offline OTs are taken with random
    // choices and corrected online.
    class PermuteShareSender : public oc::TimerAdapter
    {
        oc::u64 mN = 0, mLogN = 0, mNumColumns = 0, mNumSwitches = 0;
        Benes mBenes;

        // Two bits per switch, distilled from the OT message at generation
        // time: the message's low bit and the low bit of its AES image. Keeps
        // the offline material at 2 bits/switch instead of a 16-byte block.
        std::vector<oc::u8> mRotBits;
        oc::BitVector mRotChoices;

    public:
        // Sizes only, so the offline phase can run before the permutation is
        // known (it is derived from the sender's set).
        void init(oc::u64 n);

        // Offline: the switch OTs. Depends only on n.
        Proto setup(oc::SilentOtExtReceiver& ot, oc::PRNG& prng, Socket& chl);

        // Routes the Benes network for this permutation. Online.
        void setPermutation(std::vector<int> perm);

        // Online: outputs this party's share of pi(x).
        Proto run(oc::BitVector& share, Socket& chl);

        std::vector<int> getPerm() { return mBenes.getPerm(); };
        std::vector<int> getInvPerm() { return mBenes.getInvPerm(); };
    };

    // Holds the input vector. It is the OT *sender*.
    class PermuteShareReceiver : public oc::TimerAdapter
    {
        oc::u64 mN = 0, mLogN = 0, mNumColumns = 0, mNumSwitches = 0;

        // Four bits per switch: for each of the two OT messages, the pair of
        // bits the two wires of that switch consume.
        std::vector<std::array<std::array<oc::u8, 2>, 2>> mSotMsgs;

        void prepareCorrection(oc::u64 depth, oc::u64 permIdx,
            oc::BitVector& src,
            std::vector<std::array<std::array<oc::u8, 2>, 2>>& otMsgs,
            std::vector<std::array<oc::u8, 2>>& corrections);

    public:
        void init(oc::u64 n);

        // Offline: the switch OTs. Depends only on n.
        Proto setup(oc::SilentOtExtSender& ot, oc::PRNG& prng, Socket& chl);

        // Online: outputs this party's share of pi(inputs).
        Proto run(const oc::BitVector& inputs, oc::BitVector& outputs,
            oc::PRNG& prng, Socket& chl);
    };
}
