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

    // Correlated permutation, the whole thing precomputed.
    //
    // The permutation is sampled at random in the offline phase and the entire
    // permute + share is run there on a random vector r held by the receiver:
    // the receiver keeps a, the sender b, with a ^ b = rho(r). Online the
    // receiver only sends d = x ^ r and the sender applies rho locally, since
    // rho(d) ^ b = rho(x) ^ a. So the switch OTs, the Benes routing and both
    // network evaluations all leave the online phase.
    //
    // A caller that needs a *specific* permutation cannot use this. The PSU
    // can: all it needs is for the layout position of an item to be hidden, so
    // the sender looks up where rho sent each of its slots and runs the final
    // OT there -- from the receiver's side those positions are a uniformly
    // random subset.

    // Holds the permutation. It is the OT *receiver*: the choice bit of switch
    // i is that switch's setting, so the offline OTs are taken with random
    // choices and corrected online.
    class PermuteShareSender : public oc::TimerAdapter
    {
        oc::u64 mN = 0, mLogN = 0, mNumColumns = 0, mNumSwitches = 0;
        Benes mBenes;

        // One bit per switch, distilled from the OT message at generation
        // time, instead of a 16-byte block.
        std::vector<oc::u8> mRotBits;
        oc::BitVector mRotChoices;
        oc::BitVector mCorrShare;   // share of rho(r) from correlate()

    public:
        // Sizes only, so the offline phase can run before the permutation is
        // known (it is derived from the sender's set).
        void init(oc::u64 n);

        // Offline: the switch OTs. Depends only on n.
        Proto setup(oc::SilentOtExtReceiver& ot, oc::PRNG& prng, Socket& chl);

        // Routes the Benes network for this permutation.
        void setPermutation(std::vector<int> perm);

        // Outputs this party's share of pi(x), where x is the receiver's
        // input to the matching call.
        Proto run(oc::BitVector& share, Socket& chl);

        // Offline: samples a random permutation, runs the protocol against the
        // receiver's random vector and keeps the resulting share. Leaves the
        // permutation available through getPerm().
        Proto correlate(oc::SilentOtExtReceiver& ot, oc::PRNG& prng,
            Socket& chl);

        // Online: the receiver's blinded input d = x ^ r maps to
        // rho(d) ^ (offline share) = rho(x) ^ (receiver's offline share).
        Proto apply(oc::BitVector& share, Socket& chl);

        std::vector<int> getPerm() { return mBenes.getPerm(); };
        const std::vector<int>& getPermRef() { return mBenes.getPermRef(); };
        // Switch 0 of column 0; structurally zero (see the Waksman note).
        oc::u8 firstSwitch() { return mBenes.getSwitchesAsBitVec()[0]; };
        std::vector<int> getInvPerm() { return mBenes.getInvPerm(); };
        const std::vector<int>& getInvPermRef() { return mBenes.getInvPermRef(); };
    };

    // Holds the input vector. It is the OT *sender*.
    class PermuteShareReceiver : public oc::TimerAdapter
    {
        oc::u64 mN = 0, mLogN = 0, mNumColumns = 0, mNumSwitches = 0;

        // One bit per OT message, two messages per switch.
        std::vector<std::array<oc::u8, 2>> mSotMsgs;
        std::vector<oc::BitVector> mScratchBottom, mScratchTop;
        oc::BitVector mCorrR;       // the random vector r
        oc::BitVector mCorrShare;   // share of rho(r) from correlate()

        void prepareCorrection(oc::u64 depth, oc::u64 permIdx,
            oc::BitVector& src,
            std::vector<std::array<oc::u8, 2>>& otMsgs,
            std::vector<oc::u8>& corrections);

    public:
        void init(oc::u64 n);

        // Offline: the switch OTs. Depends only on n.
        Proto setup(oc::SilentOtExtSender& ot, oc::PRNG& prng, Socket& chl);

        // Outputs this party's share of pi(inputs).
        Proto run(const oc::BitVector& inputs, oc::BitVector& outputs,
            oc::PRNG& prng, Socket& chl);

        // Offline half of the correlation: picks the random vector r and keeps
        // the share of rho(r).
        Proto correlate(oc::SilentOtExtSender& ot, oc::PRNG& prng, Socket& chl);

        // Online: sends x ^ r and returns the offline share, which is this
        // party's share of rho(x).
        Proto apply(const oc::BitVector& inputs, oc::BitVector& outputs,
            Socket& chl);
    };
}
