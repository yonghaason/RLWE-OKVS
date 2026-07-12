#include "cpsi.h"
#include "oprf.h"
#include "band_okvs.h"
#include "../GMW/Gmw.h"

#include "cryptoTools/Common/CuckooIndex.h"
#include "cryptoTools/Common/Matrix.h"

#include <array>
#include <cstring>

using namespace std;
using namespace oc;

namespace rlweOkvs
{
    namespace
    {
        // Cuckoo hashing: 3 hash functions, no stash (as in vole-psi's CPSI).
        constexpr u64 kNumHashes = 3;

        // Band-OKVS parameters, shared with band_okvs/oprf.cpp.
        constexpr int kBandLength = 196;
        constexpr double kOkvsExpansion = 1.1;

        // Per-hash-function OPPRF key derivation, exactly as in vole-psi
        // (RsCpsi.cpp): the OPPRF key of item v under hash function j is
        // AES_{k_j ^ cuckooSeed}(v). Distinct j give distinct keys, so one
        // item inserted under all 3 hash functions yields 3 distinct
        // OKVS keys even when two of its bins collide.
        array<AES, kNumHashes> opprfHashers(block cuckooSeed)
        {
            return { AES(block(3242, 23423) ^ cuckooSeed),
                     AES(block(4534, 45654) ^ cuckooSeed),
                     AES(block(5677, 67867) ^ cuckooSeed) };
        }

        // The ss-equality compares keyBitLength = ssp + log2(numBins) bits:
        // a non-matching (pseudorandom) OPPRF output then collides with the
        // sender's bin key with probability 2^-keyBitLength per bin,
        // i.e. 2^-ssp over all numBins bins.
        u64 equalityBitLength(u64 ssp, u64 numBins)
        {
            return ssp + log2ceil(numBins);
        }

        // Length (in field elements) of the OPPRF hint OKVS. Determined by
        // the sender's 3*|Y| programmed points; both parties must compute the
        // same value.
        u64 opprfHintLength(u64 senderSize)
        {
            return static_cast<u64>(kOkvsExpansion * kNumHashes * senderSize);
        }
    }

    Proto CpsiSender::run(const vector<block>& Y, CpsiSharing& ret, Socket& chl)
    {
        if (Y.size() != mSenderSize)
            throw RTE_LOC;

        setTimePoint("CpsiSender::begin");

        // ---------------------------------------------------------------
        // Hashing: place every y under all 3 cuckoo hash functions.
        // ---------------------------------------------------------------
        block cuckooSeed;
        co_await chl.recv(cuckooSeed);

        // The receiver initialized its cuckoo table with the same
        // parameters, so numBins agrees between the parties.
        auto params = CuckooIndex<>::selectParams(mRecverSize, mSsp, 0, kNumHashes);
        u64 numBins = params.numBins();

        // Bin locations of every item under each hash function, computed with
        // the CuckooIndex machinery so they match the receiver's table
        // exactly (cf. vole-psi's SimpleIndex, which sets up mMods the same way).
        CuckooIndex<> cuckoo;
        cuckoo.mMods.resize(kNumHashes);
        for (u64 j = 0; j < kNumHashes; ++j)
            cuckoo.mMods[j] = Mod(numBins - j);

        vector<block> itemDigests(mSenderSize);
        AES(cuckooSeed).hashBlocks(Y, itemDigests);

        Matrix<u32> locations(mSenderSize, kNumHashes);
        cuckoo.computeLocations(itemDigests, locations);

        setTimePoint("CpsiSender::simple hash");

        // ---------------------------------------------------------------
        // OPPRF points: under hash function j, item y is programmed as
        //   AES_j(y)  ->  r[bin_j(y)],
        // where r[i] is this party's uniform equality key for bin i.
        // 3*|Y| points in total; the receiver's decode at a real bin
        // (holding x under hash j) queries AES_j(x) and hence recovers
        // r[bin] iff x is one of our items.
        // ---------------------------------------------------------------
        u64 keyBitLength = equalityBitLength(mSsp, numBins);
        u64 keyByteLength = divCeil(keyBitLength, 8);

        Matrix<u8> binKeys(numBins, keyByteLength, AllocType::Uninitialized);
        mPrng.get<u8>(binKeys);

        u64 numPoints = kNumHashes * mSenderSize;
        vector<block> opprfKeys(numPoints);
        vector<block> opprfValues(numPoints);
        {
            auto hashers = opprfHashers(cuckooSeed);
            for (u64 b = 0, k = 0; b < mSenderSize; ++b)
            {
                for (u64 j = 0; j < kNumHashes; ++j, ++k)
                {
                    opprfKeys[k] = hashers[j].hashBlock(Y[b]);
                    opprfValues[k] = ZeroBlock;
                    memcpy(&opprfValues[k], &binKeys(locations(b, j), 0), keyByteLength);
                }
            }
        }

        // ---------------------------------------------------------------
        // OPPRF = OPRF + OKVS hint. We evaluate the OPRF at our 3*|Y|
        // programmed keys; the interactive cost (the VOLE) is sized by the
        // receiver's numBins evaluation points. The hint encodes
        //   hint[key] = value ^ F(key)
        // over the 3*|Y| points, so its size scales with |Y|, not numBins.
        // ---------------------------------------------------------------
        OprfSender oprf;
        oprf.init(numPoints, numBins, mPrng.get());
        vector<block> F;
        co_await oprf.run(opprfKeys, F, chl);

        setTimePoint("CpsiSender::oprf");

        for (u64 k = 0; k < numPoints; ++k)
            opprfValues[k] ^= F[k];

        // Note: BandOkvs is hardwired to 128-bit values (SIMD decode paths),
        // so each hint entry ships 16 bytes even though keyByteLength (~8)
        // would suffice; vole-psi's byte-width Baxos hint is narrower per
        // entry. A 64-bit-value band OKVS would halve this term.
        band_okvs::BandOkvs okvs;
        okvs.Init(numPoints, opprfHintLength(mSenderSize), kBandLength, mPrng.get());
        vector<block> hint(okvs.Size());
        if (!okvs.Encode(opprfKeys.data(), opprfValues.data(), hint.data()))
            throw std::runtime_error(
                "OPPRF hint encoding failed: duplicate items in Y, "
                "or a band-OKVS failure event. " LOCATION);

        co_await chl.send(std::move(hint));

        setTimePoint("CpsiSender::opprf hint");

        // ---------------------------------------------------------------
        // ss-equality (GMW): share of 1(receiver's OPPRF output == r[i])
        // for every bin i, i.e. of the membership flag.
        // ---------------------------------------------------------------
        volePSI::Gmw gmw;
        if (mTimer)
            gmw.setTimer(getTimer());
        auto cir = volePSI::isZeroCircuit(keyBitLength);
        gmw.init(numBins, cir, mNumThreads, mOTeBatchSize, 1, mPrng.get());
        gmw.setInput(0, binKeys);
        co_await gmw.run(chl);

        auto shares = gmw.getOutputView(0);
        ret.mNumBins = numBins;
        ret.mFlagShares.resize(numBins);
        copy(shares.begin(), shares.begin() + ret.mFlagShares.sizeBytes(),
             ret.mFlagShares.data());

        setTimePoint("CpsiSender::gmw");
    }

    Proto CpsiReceiver::run(const vector<block>& X, CpsiSharing& ret, Socket& chl)
    {
        if (X.size() != mRecverSize)
            throw RTE_LOC;

        setTimePoint("CpsiReceiver::begin");

        // ---------------------------------------------------------------
        // Cuckoo-hash X: every item lands in one bin under one of the 3
        // hash functions; the remaining bins stay empty.
        // ---------------------------------------------------------------
        block cuckooSeed = mPrng.get();
        co_await chl.send(block(cuckooSeed));

        CuckooIndex<> cuckoo;
        cuckoo.init(mRecverSize, mSsp, 0, kNumHashes);
        // CuckooIndex::insert only reads the items; the const_cast bridges
        // its non-const span parameter.
        cuckoo.insert(oc::span<block>(const_cast<block*>(X.data()), X.size()), cuckooSeed);

        u64 numBins = cuckoo.mBins.size();

        setTimePoint("CpsiReceiver::cuckoo");

        // ---------------------------------------------------------------
        // OPPRF queries: bin i holding item x under hash function j queries
        // AES_j(x); an empty bin queries the (distinct) constant block(i, 0),
        // which the sender programmed nowhere, so its output is
        // pseudorandom and the equality below shares a 0.
        // ---------------------------------------------------------------
        vector<block> queries(numBins);
        ret.mMapping.assign(mRecverSize, ~u64(0));
        {
            auto hashers = opprfHashers(cuckooSeed);
            for (u64 i = 0; i < numBins; ++i)
            {
                auto& bin = cuckoo.mBins[i];
                if (!bin.isEmpty())
                {
                    auto b = bin.idx();
                    queries[i] = hashers[bin.hashIdx()].hashBlock(X[b]);
                    ret.mMapping[b] = i;
                }
                else
                {
                    queries[i] = block(i, 0);
                }
            }
        }
        for (u64 b = 0; b < mRecverSize; ++b)
            if (ret.mMapping[b] == ~u64(0))
                throw std::runtime_error("cuckoo insertion failed " LOCATION);

        // ---------------------------------------------------------------
        // OPPRF evaluation: OPRF at our numBins bin keys, plus the local
        // decode of the sender's hint at the same keys.
        // ---------------------------------------------------------------
        OprfReceiver oprf;
        oprf.init(numBins, kNumHashes * mSenderSize, mPrng.get());
        vector<block> F;
        co_await oprf.run(queries, F, chl);

        setTimePoint("CpsiReceiver::oprf");

        vector<block> hint;
        co_await chl.recvResize(hint);

        band_okvs::BandOkvs okvs;
        okvs.Init(numBins, opprfHintLength(mSenderSize), kBandLength);
        if (hint.size() != static_cast<u64>(okvs.Size()))
            throw std::runtime_error("unexpected OPPRF hint size " LOCATION);

        vector<block> outputs(numBins);
        okvs.Decode(queries.data(), hint.data(), outputs.data());
        for (u64 i = 0; i < numBins; ++i)
            outputs[i] ^= F[i];

        setTimePoint("CpsiReceiver::opprf decode");

        // ---------------------------------------------------------------
        // ss-equality (GMW), on the low keyBitLength bits of the OPPRF
        // outputs, against the sender's bin keys.
        // ---------------------------------------------------------------
        u64 keyBitLength = equalityBitLength(mSsp, numBins);
        u64 keyByteLength = divCeil(keyBitLength, 8);

        Matrix<u8> gmwIn(numBins, keyByteLength, AllocType::Uninitialized);
        for (u64 i = 0; i < numBins; ++i)
            memcpy(&gmwIn(i, 0), &outputs[i], keyByteLength);

        volePSI::Gmw gmw;
        if (mTimer)
            gmw.setTimer(getTimer());
        auto cir = volePSI::isZeroCircuit(keyBitLength);
        gmw.init(numBins, cir, mNumThreads, mOTeBatchSize, 0, mPrng.get());
        gmw.setInput(0, gmwIn);
        co_await gmw.run(chl);

        auto shares = gmw.getOutputView(0);
        ret.mNumBins = numBins;
        ret.mFlagShares.resize(numBins);
        copy(shares.begin(), shares.begin() + ret.mFlagShares.sizeBytes(),
             ret.mFlagShares.data());

        setTimePoint("CpsiReceiver::gmw");
    }
}
