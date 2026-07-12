#include "cpsi.h"
#include "oprf.h"
#include "band_okvs.h"
#include "../GMW/Gmw.h"

#include "macoro/when_all.h"

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
        // receiver's |X| evaluation points (its real bins -- empty bins
        // never touch the OPPRF). The hint encodes
        //   hint[key] = value ^ F(key)
        // over the 3*|Y| points, so its size scales with |Y|, not numBins.
        // ---------------------------------------------------------------
        // ---------------------------------------------------------------
        // ss-equality (GMW): share of 1(receiver's OPPRF output == r[i]) for
        // every bin i. numBins is known already, so the (input-independent)
        // triple generation is started up front on a forked channel and runs
        // concurrently with the OPRF+OPPRF phase on the base channel; only the
        // online GMW waits for the equality inputs. (Needs >=2 threads/party
        // to parallelize; with 1 thread the two just interleave.)
        // ---------------------------------------------------------------
        volePSI::Gmw gmw;
        if (mTimer)
            gmw.setTimer(getTimer());
        auto cir = volePSI::isZeroCircuit(keyBitLength);
        gmw.init(numBins, cir, mNumThreads, mOTeBatchSize, 1, mPrng.get());

        auto chlGmw = chl.fork();
        // OPPRF phase (base channel). Bound to a named local before awaiting so
        // the coroutine's by-reference captures stay alive across the join.
        auto opprfBody = [&, this]() -> Proto {
            OprfSender oprf;
            oprf.init(numPoints, mRecverSize, mPrng.get());
            vector<block> F;
            co_await oprf.run(opprfKeys, F, chl);
            setTimePoint("CpsiSender::oprf");

            for (u64 k = 0; k < numPoints; ++k)
                opprfValues[k] ^= F[k];

            // Note: BandOkvs is hardwired to 128-bit values (SIMD decode
            // paths), so each hint entry ships 16 bytes even though
            // keyByteLength (~8) would suffice; vole-psi's byte-width Baxos
            // hint is narrower. A 64-bit-value band OKVS would halve this term.
            band_okvs::BandOkvs okvs;
            okvs.Init(numPoints, opprfHintLength(mSenderSize), kBandLength,
                      mPrng.get());
            vector<block> hint(okvs.Size());
            if (!okvs.Encode(opprfKeys.data(), opprfValues.data(), hint.data()))
                throw std::runtime_error(
                    "OPPRF hint encoding failed: duplicate items in Y, "
                    "or a band-OKVS failure event. " LOCATION);
            co_await chl.send(std::move(hint));
            setTimePoint("CpsiSender::opprf hint");
        };
        auto opprf = opprfBody();
        auto tg = gmw.generateTriple(chlGmw);
        auto both = co_await macoro::when_all_ready(std::move(opprf), std::move(tg));
        std::get<0>(both).result();
        std::get<1>(both).result();

        gmw.setInput(0, binKeys);
        co_await gmw.run(chl);  // triples ready -> online only

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
        // OPPRF queries: only the |X| REAL bins touch the OPPRF -- the bin
        // holding item x under hash function j queries AES_j(x). Empty bins
        // need no OPPRF output at all: they enter the equality below with
        // private uniform randomness (and hence share a 0 up to a
        // 2^-keyBitLength collision, the same bound as a non-member).
        // vole-psi instead evaluates at all numBins bins, dummies included;
        // the interactive OPRF here is sized by |X|, not numBins.
        // ---------------------------------------------------------------
        vector<block> queries(mRecverSize);
        vector<u64> queryBin(mRecverSize);
        ret.mMapping.assign(mRecverSize, ~u64(0));
        {
            auto hashers = opprfHashers(cuckooSeed);
            u64 k = 0;
            for (u64 i = 0; i < numBins; ++i)
            {
                auto& bin = cuckoo.mBins[i];
                if (!bin.isEmpty())
                {
                    auto b = bin.idx();
                    queries[k] = hashers[bin.hashIdx()].hashBlock(X[b]);
                    queryBin[k] = i;
                    ret.mMapping[b] = i;
                    ++k;
                }
            }
            if (k != mRecverSize)
                throw std::runtime_error("cuckoo insertion failed " LOCATION);
        }

        // ---------------------------------------------------------------
        // ss-equality (GMW), on the low keyBitLength bits of the OPPRF
        // outputs, against the sender's bin keys. The equality covers every
        // bin -- empty ones with fresh randomness -- so the sender learns
        // nothing about which bins are occupied. numBins is known, so the
        // (input-independent) triple generation runs on a forked channel
        // concurrently with the OPRF+OPPRF phase on the base channel; only
        // the online GMW waits for the outputs.
        // ---------------------------------------------------------------
        u64 keyBitLength = equalityBitLength(mSsp, numBins);
        u64 keyByteLength = divCeil(keyBitLength, 8);

        volePSI::Gmw gmw;
        if (mTimer)
            gmw.setTimer(getTimer());
        auto cir = volePSI::isZeroCircuit(keyBitLength);
        gmw.init(numBins, cir, mNumThreads, mOTeBatchSize, 0, mPrng.get());

        // OPPRF evaluation (base channel): OPRF at our |X| bin keys plus the
        // local decode of the sender's hint. Fills `outputs`.
        vector<block> outputs(mRecverSize);
        auto chlGmw = chl.fork();
        auto opprfBody = [&, this]() -> Proto {
            OprfReceiver oprf;
            oprf.init(mRecverSize, kNumHashes * mSenderSize, mPrng.get());
            vector<block> F;
            co_await oprf.run(queries, F, chl);
            setTimePoint("CpsiReceiver::oprf");

            vector<block> hint;
            co_await chl.recvResize(hint);

            band_okvs::BandOkvs okvs;
            okvs.Init(mRecverSize, opprfHintLength(mSenderSize), kBandLength);
            if (hint.size() != static_cast<u64>(okvs.Size()))
                throw std::runtime_error("unexpected OPPRF hint size " LOCATION);

            okvs.Decode(queries.data(), hint.data(), outputs.data());
            for (u64 k = 0; k < mRecverSize; ++k)
                outputs[k] ^= F[k];
            setTimePoint("CpsiReceiver::opprf decode");
        };
        auto opprf = opprfBody();
        auto tg = gmw.generateTriple(chlGmw);
        auto both = co_await macoro::when_all_ready(std::move(opprf), std::move(tg));
        std::get<0>(both).result();
        std::get<1>(both).result();

        Matrix<u8> gmwIn(numBins, keyByteLength, AllocType::Uninitialized);
        mPrng.get<u8>(gmwIn);
        for (u64 k = 0; k < mRecverSize; ++k)
            memcpy(&gmwIn(queryBin[k], 0), &outputs[k], keyByteLength);

        gmw.setInput(0, gmwIn);
        co_await gmw.run(chl);  // triples ready -> online only

        auto shares = gmw.getOutputView(0);
        ret.mNumBins = numBins;
        ret.mFlagShares.resize(numBins);
        copy(shares.begin(), shares.begin() + ret.mFlagShares.sizeBytes(),
             ret.mFlagShares.data());

        setTimePoint("CpsiReceiver::gmw");
    }
}
