#pragma once
#include "cryptoTools/Common/Defines.h"
#include "cryptoTools/Common/block.h"
#include "cryptoTools/Common/BitVector.h"
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Common/Timer.h"
#include "coproto/coproto.h"

#include <vector>

namespace rlweOkvs
{
    using Proto = coproto::task<>;
    using Socket = coproto::Socket;

    // Circuit-PSI-based secret-shared private membership test (ss-PMT).
    //
    // This is the membership-share core of vole-psi's RsCpsi (the PSTY19-style
    // circuit PSI; see volePSI/RsCpsi.{h,cpp}), with two deliberate changes:
    //
    //  1. Every OKVS is the band OKVS of this repository:
    //     - the OPRF is the VOLE OPRF of band_okvs/oprf.{h,cpp} (RR22 style,
    //       band OKVS instead of the RR22 Baxos), and
    //     - the OPPRF hint is a band OKVS over the sender's 3*|Y| programmed
    //       points (vole-psi likewise programs exactly 3*n_s points -- the
    //       hint never scales with the number of cuckoo bins).
    //  2. The receiver runs the OPPRF only at its |X| REAL bins, not at all
    //     numBins bins as vole-psi does: an empty bin needs no OPPRF output,
    //     it can feed private uniform randomness straight into the equality
    //     (it then shares a 0 up to the same 2^-ssp collision bound as a
    //     non-member). The interactive cost therefore scales with |X| and
    //     3*|Y| only; nothing interactive scales with numBins except the
    //     equality itself, which must cover every bin so that the sender
    //     learns nothing about which bins are occupied.
    //  3. The payload-sharing machinery of RsCpsi (associated values and their
    //     Gmw muxing) is omitted: the parties obtain XOR shares of the
    //     membership flags only, which is all an ss-PMT needs.
    //
    // Roles (identical to vole-psi, so the two are easy to diff):
    //  - the RECEIVER holds X and cuckoo-hashes it into numBins bins
    //    (3 hash functions, no stash, numBins ~ 1.27|X|);
    //  - the SENDER holds Y and inserts every item under all 3 hash
    //    functions, i.e. programs 3*|Y| OPPRF points.
    //
    // Output: for every cuckoo bin i the parties hold bits s_i (sender) and
    // r_i (receiver) with s_i ^ r_i = 1(the receiver's item in bin i is in Y).
    // The receiver knows which bin holds which of its items (mMapping);
    // an empty bin shares a 0, and a real non-member bin shares a 0 up to the
    // 2^-ssp equality-collision probability.

    struct CpsiSharing
    {
        // One share bit per cuckoo bin.
        oc::BitVector mFlagShares;

        // Receiver only: mMapping[b] is the bin holding X[b]
        // (always valid -- cuckoo insertion failure throws).
        std::vector<oc::u64> mMapping;

        oc::u64 mNumBins = 0;
    };

    class CpsiSender : public oc::TimerAdapter
    {
        oc::u64 mSenderSize = 0;
        oc::u64 mRecverSize = 0;
        oc::u64 mSsp = 40;
        oc::u64 mNumThreads = 1;
        oc::u64 mOTeBatchSize = 1ull << 19;
        oc::PRNG mPrng;

    public:
        void init(oc::u64 senderSize, oc::u64 recverSize, oc::block seed,
                  oc::u64 statSecParam = 40, oc::u64 numThreads = 1)
        {
            mSenderSize = senderSize;
            mRecverSize = recverSize;
            mSsp = statSecParam;
            mNumThreads = numThreads;
            mPrng.SetSeed(seed);
        }

        // Y: the sender's set. Items must be distinct (they key an OKVS).
        Proto run(const std::vector<oc::block>& Y, CpsiSharing& out, Socket& chl);
    };

    class CpsiReceiver : public oc::TimerAdapter
    {
        oc::u64 mSenderSize = 0;
        oc::u64 mRecverSize = 0;
        oc::u64 mSsp = 40;
        oc::u64 mNumThreads = 1;
        oc::u64 mOTeBatchSize = 1ull << 19;
        oc::PRNG mPrng;

    public:
        void init(oc::u64 recverSize, oc::u64 senderSize, oc::block seed,
                  oc::u64 statSecParam = 40, oc::u64 numThreads = 1)
        {
            mRecverSize = recverSize;
            mSenderSize = senderSize;
            mSsp = statSecParam;
            mNumThreads = numThreads;
            mPrng.SetSeed(seed);
        }

        // X: the receiver's set. Items must be distinct (they are cuckoo-hashed).
        Proto run(const std::vector<oc::block>& X, CpsiSharing& out, Socket& chl);
    };
}
