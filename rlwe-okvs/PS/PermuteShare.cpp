#include "PermuteShare.h"

#include <cmath>
#include <cstring>
#include <numeric>


using namespace oc;

namespace rlweOkvs {
namespace {
// A switch consumes a single bit from its OT message.
inline u8 msgBit(const block &m) { return (u8)(m.mData[0] & 1); }

// Silent OT in chunks, so the block buffers never hold all mNumSwitches at
// once; each chunk is distilled to bits and dropped.
constexpr u64 kOtChunk = 1ull << 22;

// BitVector packs bits LSB-first. Walking a byte at a time and shifting keeps
// the per-bit work to a shift and a mask, instead of the index arithmetic and
// bounds check operator[] does -- these loops run once per switch, tens of
// millions of times.
// The bits are random, so every branch on one is a coin flip the predictor
// cannot win; these loops stay branchless.
inline u8 bitAt(const u8 *bytes, u64 i) { return (bytes[i >> 3] >> (i & 7)) & 1; }
inline void orBit(u8 *bytes, u64 i, u8 v) {
  bytes[i >> 3] |= (u8)((v & 1) << (i & 7));  // buffer starts zeroed
}
}  // namespace

u64 permuteShareSwitchCount(u64 n) {
  if (n < 2) return 0;
  const u64 logN = (u64)std::ceil(std::log2((double)n));
  return (2 * logN - 1) * (n / 2);
}

// ---------------------------------------------------------------- sender ---

void PermuteShareSender::init(u64 n) {
  mN = n;
  mLogN = (u64)std::ceil(std::log2((double)n));
  mNumColumns = 2 * mLogN - 1;
  mNumSwitches = mNumColumns * (mN / 2);
}

void PermuteShareSender::setPermutation(std::vector<int> perm) {
  mBenes.initialize((int)mN, std::move(perm));
  mBenes.genBenesRoute();
  setTimePoint("PermuteShareSender::Benes routing");
}

Proto PermuteShareSender::setup(SilentOtExtReceiver &ot, PRNG &prng,
                                Socket &chl) {
  mRotBits.resize(mNumSwitches);
  mRotChoices.resize(mNumSwitches);

  std::vector<block> msgs;
  BitVector choices;
  for (u64 begin = 0; begin < mNumSwitches; begin += kOtChunk) {
    const u64 len = std::min<u64>(kOtChunk, mNumSwitches - begin);
    msgs.resize(len);
    choices.resize(len);
    co_await ot.silentReceive(choices, msgs, prng, chl, oc::OTType::Random);
    for (u64 i = 0; i < len; ++i) {
      mRotBits[begin + i] = msgBit(msgs[i]);
      mRotChoices[begin + i] = choices[i];
    }
  }
  setTimePoint("PermuteShareSender::Setup (switch OTs)");
}

Proto PermuteShareSender::run(BitVector &share, Socket &chl) {
  // Steer each OT onto this switch's actual setting.
  setTimePoint("PS.S: online begin");
  BitVector switches = mBenes.getSwitchesAsBitVec();
  setTimePoint("PS.S: switches as bitvec");
  BitVector bitCorrection = switches;
  bitCorrection ^= mRotChoices;
  co_await chl.send(bitCorrection);

  BitVector recvCorr(mNumSwitches);
  co_await chl.recv(recvCorr);
  setTimePoint("PS.S: recv corrections");

  std::vector<u8> recvMsg(mNumSwitches);
  {
    const u8 *sw = switches.data();
    const u8 *cc = recvCorr.data();
    for (u64 i = 0; i < mNumSwitches; ++i) {
      const u8 m = (u8)(0 - bitAt(sw, i));  // 0x00 or 0xff
      recvMsg[i] = (u8)(mRotBits[i] ^ (bitAt(cc, i) & m));
    }
  }

  setTimePoint("PS.S: apply corrections");

  share.resize(mN);
  co_await chl.recv(share);

  // recvMsg is already laid out column * (mN/2) + switch, which is what the
  // evaluation indexes -- no matrix to build.
  mBenes.benesMaskedEval(share, recvMsg);
  setTimePoint("PS.S: masked eval");
}


Proto PermuteShareSender::correlate(SilentOtExtReceiver &ot, PRNG &prng,
                                    Socket &chl) {
  co_await setup(ot, prng, chl);

  // A random permutation is as good as any: the caller only needs positions
  // hidden, and it can look up where rho sent each of them afterwards.
  std::vector<int> perm(mN);
  std::iota(perm.begin(), perm.end(), 0);
  for (u64 i = mN - 1; i > 0; --i) {
    std::swap(perm[i], perm[prng.get<u64>() % (i + 1)]);
  }
  setPermutation(std::move(perm));

  co_await run(mCorrShare, chl);
  setTimePoint("PermuteShareSender::Correlate");
}

Proto PermuteShareSender::apply(BitVector &share, Socket &chl) {
  BitVector d(mN);
  co_await chl.recv(d);

  // rho(d) ^ b = rho(x) ^ rho(r) ^ b = rho(x) ^ a.
  const auto &perm = mBenes.getPermRef();
  share.resize(mN);
  for (u64 i = 0; i < mN; ++i) {
    share[i] = d[perm[i]] ^ mCorrShare[i];
  }
  setTimePoint("PermuteShareSender::Apply");
}

// -------------------------------------------------------------- receiver ---

void PermuteShareReceiver::init(u64 n) {
  mN = n;
  mLogN = (u64)std::ceil(std::log2((double)n));
  mNumColumns = 2 * mLogN - 1;
  mNumSwitches = mNumColumns * (mN / 2);
}

Proto PermuteShareReceiver::setup(SilentOtExtSender &ot, PRNG &prng,
                                  Socket &chl) {
  mSotMsgs.resize(mNumSwitches);

  std::vector<std::array<block, 2>> msgs;
  for (u64 begin = 0; begin < mNumSwitches; begin += kOtChunk) {
    const u64 len = std::min<u64>(kOtChunk, mNumSwitches - begin);
    msgs.resize(len);
    co_await ot.silentSend(msgs, prng, chl);
    for (u64 i = 0; i < len; ++i) {
      mSotMsgs[begin + i][0] = msgBit(msgs[i][0]);
      mSotMsgs[begin + i][1] = msgBit(msgs[i][1]);
    }
  }
  setTimePoint("PermuteShareReceiver::Setup (switch OTs)");
}

Proto PermuteShareReceiver::run(const BitVector &inputs, BitVector &outputs,
                                PRNG &prng, Socket &chl) {
  BitVector masks(mN);
  masks.randomize(prng);
  outputs = masks;  // the input-side masks; permuted in place below

  setTimePoint("PS.R: online begin");
  BitVector bitCorrection(mNumSwitches);
  co_await chl.recv(bitCorrection);
  setTimePoint("PS.R: recv bitcorrection");

  // Swap the message pair wherever the sender's switch disagreed with the
  // random choice bit of the offline OT. Done in place: the offline material
  // is consumed by this run.
  auto &sot = mSotMsgs;
  {
    const u8 *bc = bitCorrection.data();
    for (u64 i = 0; i < mNumSwitches; ++i) {
      const u8 m = (u8)(0 - bitAt(bc, i));
      const u8 t = (u8)((sot[i][0] ^ sot[i][1]) & m);
      sot[i][0] ^= t;
      sot[i][1] ^= t;
    }
  }

  setTimePoint("PS.R: copy+swap sot");

  // Depth-first recursion: one scratch buffer per level suffices, and sizing
  // them up front removes the two allocations every call was making.
  mScratchBottom.resize(mLogN + 2);
  mScratchTop.resize(mLogN + 2);
  for (u64 i = 0; i < mScratchBottom.size(); ++i) {
    mScratchBottom[i].reserve(mN);
    mScratchTop[i].reserve(mN);
  }

  std::vector<u8> corrections(mNumSwitches);
  BitVector permuted = masks;
  prepareCorrection(0, 0, permuted, sot, corrections);
  setTimePoint("PS.R: prepare corrections");

  BitVector corrBits(mNumSwitches);
  {
    u8 *b = corrBits.data();
    for (u64 i = 0; i < mNumSwitches; ++i) orBit(b, i, corrections[i]);
  }
  setTimePoint("PS.R: pack corrections");
  co_await chl.send(corrBits);

  // The sender evaluates the network on (mask ^ input), so its output and the
  // permuted masks are shares of pi(input).
  BitVector blinded = masks;
  blinded ^= inputs;
  co_await chl.send(blinded);

  outputs = permuted;
  setTimePoint("PermuteShareReceiver::Online");
}

  void PermuteShareReceiver::prepareCorrection(
    u64 depth, u64 permIdx, oc::BitVector &src,
    std::vector<std::array<u8, 2>> &otMsgs,
    std::vector<u8> &corrections)
  {
    // ot message M0 = m0 ^ w0 || m1 ^ w1
    //  for each switch: top wire m0 w0 - bottom wires m1, w1
    //  M1 = m0 ^ w1 || m1 ^ w0
    int coDepth = mLogN - depth;
    int levels = 2 * coDepth - 1;

    int subNetSize = src.size();

    // One-bit messages: the receiver picks the outgoing masks as
    // w0 = M0 ^ m0 and w1 = M0 ^ m1, so both wires consume the same bit and
    // the corrected M1 is M0 ^ m0 ^ m1. That is what the XOR of a switch's two
    // outputs being invariant buys -- half the message, half the correction.
    u8 m0, m1, w0, w1, M0, M1;
    int baseIdx;

    oc::BitVector &bottom1 = mScratchBottom[depth];
    oc::BitVector &top1 = mScratchTop[depth];
    bottom1.resize(0);
    top1.resize(0);

    if (subNetSize == 2)
    {
      if (coDepth == 1)
      {
        baseIdx = depth * (mN / 2) + permIdx;
        m0 = (u8)src[0];
        m1 = (u8)src[1];
        M0 = otMsgs[baseIdx][0];
      w0 = M0 ^ m0;
      w1 = M0 ^ m1;
      M1 = (u8)(M0 ^ m0 ^ m1);
      corrections[baseIdx] = (u8)(otMsgs[baseIdx][1] ^ M1);
        src[0] = w0;
        src[1] = w1;
      }
      else
      {
        baseIdx = (depth + 1) * (mN / 2) + permIdx;
        m0 = (u8)src[0];
        m1 = (u8)src[1];
        M0 = otMsgs[baseIdx][0];
      w0 = M0 ^ m0;
      w1 = M0 ^ m1;
      M1 = (u8)(M0 ^ m0 ^ m1);
      corrections[baseIdx] = (u8)(otMsgs[baseIdx][1] ^ M1);
        src[0] = w0;
        src[1] = w1;
      }
      return;
    }

    if (subNetSize == 3)
    {
      baseIdx = depth * (mN / 2) + permIdx;
      m0 = (u8)src[0];
      m1 = (u8)src[1];
      M0 = otMsgs[baseIdx][0];
      w0 = M0 ^ m0;
      w1 = M0 ^ m1;
      M1 = (u8)(M0 ^ m0 ^ m1);
      corrections[baseIdx] = (u8)(otMsgs[baseIdx][1] ^ M1);
      src[0] = w0;
      src[1] = w1;

      baseIdx = (depth + 1) * (mN / 2) + permIdx;
      m0 = (u8)src[1];
      m1 = (u8)src[2];
      M0 = otMsgs[baseIdx][0];
      w0 = M0 ^ m0;
      w1 = M0 ^ m1;
      M1 = (u8)(M0 ^ m0 ^ m1);
      corrections[baseIdx] = (u8)(otMsgs[baseIdx][1] ^ M1);
      src[1] = w0;
      src[2] = w1;

      baseIdx = (depth + 2) * (mN / 2) + permIdx;
      m0 = (u8)src[0];
      m1 = (u8)src[1];
      M0 = otMsgs[baseIdx][0];
      w0 = M0 ^ m0;
      w1 = M0 ^ m1;
      M1 = (u8)(M0 ^ m0 ^ m1);
      corrections[baseIdx] = (u8)(otMsgs[baseIdx][1] ^ M1);
      src[0] = w0;
      src[1] = w1;
      return;
    }

    // partea superioara. The receiver's masks do not swap -- it does not know
    // the switch settings -- so the split is already branchless; presizing
    // just avoids the pushBack bookkeeping.
    bottom1.resize(subNetSize / 2);
    top1.resize(subNetSize / 2 + (subNetSize % 2));
    for (int i = 0; i < subNetSize - 1; i += 2)
    {
      baseIdx = (depth) * (mN / 2) + permIdx + i / 2;
      m0 = (u8)src[i];
      m1 = (u8)src[i ^ 1];
      M0 = otMsgs[baseIdx][0];
      w0 = M0 ^ m0;
      w1 = M0 ^ m1;
      M1 = (u8)(M0 ^ m0 ^ m1);
      corrections[baseIdx] = (u8)(otMsgs[baseIdx][1] ^ M1);
      src[i] = w0;
      src[i ^ 1] = w1;

      bottom1[i / 2] = (u8)src[i];
      top1[i / 2] = (u8)src[i ^ 1];
    }

    if (subNetSize % 2 == 1)
    {
      top1[subNetSize / 2] = (u8)src[subNetSize - 1];
    }

    prepareCorrection(depth + 1, permIdx + subNetSize / 4, top1, otMsgs, corrections);
    prepareCorrection(depth + 1, permIdx, bottom1, otMsgs, corrections);

    // partea inferioara
    for (int i = 0; i < subNetSize - 1; i += 2)
    {
      baseIdx = (depth + levels - 1) * (mN / 2) + permIdx + i / 2;
      m1 = (u8)top1[i / 2];
      m0 = (u8)bottom1[i / 2];
      M0 = otMsgs[baseIdx][0];
      w0 = M0 ^ m0;
      w1 = M0 ^ m1;
      M1 = (u8)(M0 ^ m0 ^ m1);
      corrections[baseIdx] = (u8)(otMsgs[baseIdx][1] ^ M1);
      src[i] = w0;
      src[i ^ 1] = w1;
    }

    int idx = int(ceil(subNetSize * 0.5));
    if (subNetSize % 2 == 1)
    {
      src[subNetSize - 1] = top1[idx - 1];
    }
  }

Proto PermuteShareReceiver::correlate(SilentOtExtSender &ot, PRNG &prng,
                                      Socket &chl) {
  co_await setup(ot, prng, chl);

  mCorrR.resize(mN);
  mCorrR.randomize(prng);
  co_await run(mCorrR, mCorrShare, prng, chl);
  setTimePoint("PermuteShareReceiver::Correlate");
}

Proto PermuteShareReceiver::apply(const BitVector &inputs, BitVector &outputs,
                                  Socket &chl) {
  BitVector d = mCorrR;
  d ^= inputs;
  co_await chl.send(d);
  outputs = mCorrShare;
  setTimePoint("PermuteShareReceiver::Apply");
}

}  // namespace rlweOkvs