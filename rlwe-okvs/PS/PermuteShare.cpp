#include "PermuteShare.h"

#include <cmath>
#include <cstring>

#include "cryptoTools/Crypto/AES.h"

using namespace oc;

namespace rlweOkvs {
namespace {
// The two bits a switch consumes from one OT message.
inline std::array<u8, 2> msgBits(const block &m, const AES &aes) {
  return {(u8)(m.mData[0] & 1), (u8)(aes.ecbEncBlock(m).mData[0] & 1)};
}

// Silent OT in chunks, so the block buffers never hold all mNumSwitches at
// once; each chunk is distilled to bits and dropped.
constexpr u64 kOtChunk = 1ull << 22;
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
  AES aes(oc::ZeroBlock);
  mRotBits.resize(mNumSwitches * 2);
  mRotChoices.resize(mNumSwitches);

  std::vector<block> msgs;
  BitVector choices;
  for (u64 begin = 0; begin < mNumSwitches; begin += kOtChunk) {
    const u64 len = std::min<u64>(kOtChunk, mNumSwitches - begin);
    msgs.resize(len);
    choices.resize(len);
    co_await ot.silentReceive(choices, msgs, prng, chl, oc::OTType::Random);
    for (u64 i = 0; i < len; ++i) {
      const auto b = msgBits(msgs[i], aes);
      mRotBits[2 * (begin + i)] = b[0];
      mRotBits[2 * (begin + i) + 1] = b[1];
      mRotChoices[begin + i] = choices[i];
    }
  }
  setTimePoint("PermuteShareSender::Setup (switch OTs)");
}

Proto PermuteShareSender::run(BitVector &share, Socket &chl) {
  // Steer each OT onto this switch's actual setting.
  BitVector switches = mBenes.getSwitchesAsBitVec();
  BitVector bitCorrection = switches;
  bitCorrection ^= mRotChoices;
  co_await chl.send(bitCorrection);

  std::vector<BitVector> recvCorr(2);
  recvCorr[0].resize(mNumSwitches);
  recvCorr[1].resize(mNumSwitches);
  co_await chl.recv(recvCorr[0]);
  co_await chl.recv(recvCorr[1]);

  std::vector<std::array<u8, 2>> recvMsg(mNumSwitches);
  for (u64 i = 0; i < mNumSwitches; ++i) {
    u8 a = mRotBits[2 * i], b = mRotBits[2 * i + 1];
    if (switches[i]) {
      a ^= (u8)recvCorr[0][i];
      b ^= (u8)recvCorr[1][i];
    }
    recvMsg[i] = {a, b};
  }

  share.resize(mN);
  co_await chl.recv(share);

  std::vector<std::vector<std::array<u8, 2>>> matrix(mNumColumns);
  u64 ctr = 0;
  for (u64 i = 0; i < mNumColumns; ++i) {
    matrix[i].resize(mN);
    for (u64 j = 0; j < mN / 2; ++j) matrix[i][j] = recvMsg[ctr++];
  }
  mBenes.benesMaskedEval(share, matrix);
  setTimePoint("PermuteShareSender::Online");
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
  AES aes(oc::ZeroBlock);
  mSotMsgs.resize(mNumSwitches);

  std::vector<std::array<block, 2>> msgs;
  for (u64 begin = 0; begin < mNumSwitches; begin += kOtChunk) {
    const u64 len = std::min<u64>(kOtChunk, mNumSwitches - begin);
    msgs.resize(len);
    co_await ot.silentSend(msgs, prng, chl);
    for (u64 i = 0; i < len; ++i) {
      mSotMsgs[begin + i][0] = msgBits(msgs[i][0], aes);
      mSotMsgs[begin + i][1] = msgBits(msgs[i][1], aes);
    }
  }
  setTimePoint("PermuteShareReceiver::Setup (switch OTs)");
}

Proto PermuteShareReceiver::run(const BitVector &inputs, BitVector &outputs,
                                PRNG &prng, Socket &chl) {
  BitVector masks(mN);
  masks.randomize(prng);
  outputs = masks;  // the input-side masks; permuted in place below

  BitVector bitCorrection(mNumSwitches);
  co_await chl.recv(bitCorrection);

  // Swap the message pair wherever the sender's switch disagreed with the
  // random choice bit of the offline OT.
  auto sot = mSotMsgs;
  for (u64 i = 0; i < mNumSwitches; ++i) {
    if (bitCorrection[i]) std::swap(sot[i][0], sot[i][1]);
  }

  std::vector<std::array<u8, 2>> corrections(mNumSwitches);
  BitVector permuted = masks;
  prepareCorrection(0, 0, permuted, sot, corrections);

  std::vector<BitVector> corrBits(2);
  corrBits[0].resize(mNumSwitches);
  corrBits[1].resize(mNumSwitches);
  for (u64 i = 0; i < mNumSwitches; ++i) {
    corrBits[0][i] = corrections[i][0];
    corrBits[1][i] = corrections[i][1];
  }
  co_await chl.send(corrBits[0]);
  co_await chl.send(corrBits[1]);

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
    std::vector<std::array<std::array<u8, 2>, 2>> &otMsgs,
    std::vector<std::array<u8, 2>> &corrections)
  {
    // ot message M0 = m0 ^ w0 || m1 ^ w1
    //  for each switch: top wire m0 w0 - bottom wires m1, w1
    //  M1 = m0 ^ w1 || m1 ^ w0
    int coDepth = mLogN - depth;
    int levels = 2 * coDepth - 1;

    int subNetSize = src.size();

    u8 m0, m1, w0, w1, M0[2], M1[2], corrMsg[2];
    std::array<u8, 2> temp;
    int baseIdx;

    oc::BitVector bottom1;
    oc::BitVector top1;

    if (subNetSize == 2)
    {
      if (coDepth == 1)
      {
        baseIdx = depth * (mN / 2) + permIdx;
        m0 = (u8)src[0];
        m1 = (u8)src[1];
        temp = otMsgs[baseIdx][0];
        memcpy(M0, temp.data(), sizeof(M0));
        w0 = M0[0] ^ m0;
        w1 = M0[1] ^ m1;
        temp = otMsgs[baseIdx][1];
        memcpy(M1, temp.data(), sizeof(M1));
        corrMsg[0] = M1[0] ^ m0 ^ w1;
        corrMsg[1] = M1[1] ^ m1 ^ w0;
        corrections[baseIdx] = {corrMsg[0], corrMsg[1]};
        M1[0] = m0 ^ w1;
        M1[1] = m1 ^ w0;
        otMsgs[baseIdx][1] = {M1[0], M1[1]};
        src[0] = w0;
        src[1] = w1;
      }
      else
      {
        baseIdx = (depth + 1) * (mN / 2) + permIdx;
        m0 = (u8)src[0];
        m1 = (u8)src[1];
        temp = otMsgs[baseIdx][0];
        memcpy(M0, temp.data(), sizeof(M0));
        w0 = M0[0] ^ m0;
        w1 = M0[1] ^ m1;
        temp = otMsgs[baseIdx][1];
        memcpy(M1, temp.data(), sizeof(M1));
        corrMsg[0] = M1[0] ^ m0 ^ w1;
        corrMsg[1] = M1[1] ^ m1 ^ w0;
        corrections[baseIdx] = {corrMsg[0], corrMsg[1]};
        M1[0] = m0 ^ w1;
        M1[1] = m1 ^ w0;
        otMsgs[baseIdx][1] = {M1[0], M1[1]};
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
      temp = otMsgs[baseIdx][0];
      memcpy(M0, temp.data(), sizeof(M0));
      w0 = M0[0] ^ m0;
      w1 = M0[1] ^ m1;
      temp = otMsgs[baseIdx][1];
      memcpy(M1, temp.data(), sizeof(M1));
      corrMsg[0] = M1[0] ^ m0 ^ w1;
      corrMsg[1] = M1[1] ^ m1 ^ w0;
      corrections[baseIdx] = {corrMsg[0], corrMsg[1]};
      M1[0] = m0 ^ w1;
      M1[1] = m1 ^ w0;
      otMsgs[baseIdx][1] = {M1[0], M1[1]};
      src[0] = w0;
      src[1] = w1;

      baseIdx = (depth + 1) * (mN / 2) + permIdx;
      m0 = (u8)src[1];
      m1 = (u8)src[2];
      temp = otMsgs[baseIdx][0];
      memcpy(M0, temp.data(), sizeof(M0));
      w0 = M0[0] ^ m0;
      w1 = M0[1] ^ m1;
      temp = otMsgs[baseIdx][1];
      memcpy(M1, temp.data(), sizeof(M1));
      corrMsg[0] = M1[0] ^ m0 ^ w1;
      corrMsg[1] = M1[1] ^ m1 ^ w0;
      corrections[baseIdx] = {corrMsg[0], corrMsg[1]};
      M1[0] = m0 ^ w1;
      M1[1] = m1 ^ w0;
      otMsgs[baseIdx][1] = {M1[0], M1[1]};
      src[1] = w0;
      src[2] = w1;

      baseIdx = (depth + 2) * (mN / 2) + permIdx;
      m0 = (u8)src[0];
      m1 = (u8)src[1];
      temp = otMsgs[baseIdx][0];
      memcpy(M0, temp.data(), sizeof(M0));
      w0 = M0[0] ^ m0;
      w1 = M0[1] ^ m1;
      temp = otMsgs[baseIdx][1];
      memcpy(M1, temp.data(), sizeof(M1));
      corrMsg[0] = M1[0] ^ m0 ^ w1;
      corrMsg[1] = M1[1] ^ m1 ^ w0;
      corrections[baseIdx] = {corrMsg[0], corrMsg[1]};
      M1[0] = m0 ^ w1;
      M1[1] = m1 ^ w0;
      otMsgs[baseIdx][1] = {M1[0], M1[1]};
      src[0] = w0;
      src[1] = w1;
      return;
    }

    // partea superioara
    for (int i = 0; i < subNetSize - 1; i += 2)
    {
      baseIdx = (depth) * (mN / 2) + permIdx + i / 2;
      m0 = (u8)src[i];
      m1 = (u8)src[i ^ 1];
      temp = otMsgs[baseIdx][0];
      memcpy(M0, temp.data(), sizeof(M0));
      w0 = M0[0] ^ m0;
      w1 = M0[1] ^ m1;
      temp = otMsgs[baseIdx][1];
      memcpy(M1, temp.data(), sizeof(M1));
      corrMsg[0] = M1[0] ^ m0 ^ w1;
      corrMsg[1] = M1[1] ^ m1 ^ w0;
      corrections[baseIdx] = {corrMsg[0], corrMsg[1]};
      M1[0] = m0 ^ w1;
      M1[1] = m1 ^ w0;
      otMsgs[baseIdx][1] = {M1[0], M1[1]};
      src[i] = w0;
      src[i ^ 1] = w1;

      bottom1.pushBack(src[i]);
      top1.pushBack(src[i ^ 1]);
    }

    if (subNetSize % 2 == 1)
    {
      top1.pushBack(src[subNetSize - 1]);
    }

    prepareCorrection(depth + 1, permIdx + subNetSize / 4, top1, otMsgs, corrections);
    prepareCorrection(depth + 1, permIdx, bottom1, otMsgs, corrections);

    // partea inferioara
    for (int i = 0; i < subNetSize - 1; i += 2)
    {
      baseIdx = (depth + levels - 1) * (mN / 2) + permIdx + i / 2;
      m1 = (u8)top1[i / 2];
      m0 = (u8)bottom1[i / 2];
      temp = otMsgs[baseIdx][0];
      memcpy(M0, temp.data(), sizeof(M0));
      w0 = M0[0] ^ m0;
      w1 = M0[1] ^ m1;
      temp = otMsgs[baseIdx][1];
      memcpy(M1, temp.data(), sizeof(M1));
      corrMsg[0] = M1[0] ^ m0 ^ w1;
      corrMsg[1] = M1[1] ^ m1 ^ w0;
      corrections[baseIdx] = {corrMsg[0], corrMsg[1]};
      M1[0] = m0 ^ w1;
      M1[1] = m1 ^ w0;
      otMsgs[baseIdx][1] = {M1[0], M1[1]};
      src[i] = w0;
      src[i ^ 1] = w1;
    }

    int idx = int(ceil(subNetSize * 0.5));
    if (subNetSize % 2 == 1)
    {
      src[subNetSize - 1] = top1[idx - 1];
    }
  }
}  // namespace rlweOkvs
