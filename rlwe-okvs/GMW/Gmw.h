// Adapted from the GMW implementation in https://github.com/ladnir/volepsi.
#pragma once

#include "Defines.h"
#include "Circuit.h"
#include "SilentTripleGen.h"
#include <list>
#include <cryptoTools/Network/Channel.h>
#include <cryptoTools/Common/Matrix.h>

using namespace std;

namespace volePSI
{
    enum class OtExtType
    {
        IKNP,
        Silent,
        InsecureMock
    };
    class Gmw : public oc::TimerAdapter
    {
    public:

        struct Debug
        {
            bool mDebug = false;
            std::vector<block> mA, mB, mC, mD;
            std::list<std::array<std::vector<block>,2>> mU, mW;
            Matrix<int> mVals;
            Matrix<block> mWords;
        };

        Debug mO;

        BetaCircuit::LevelizeType mLevelize = BetaCircuit::LevelizeType::Reorder;

        u64 mN = 0, mNumOts = 0, mIdx, mNumThreads, mBatchSize;
        OtExtType mOtExtType;

        Matrix<block> mWords;
        u64 mRoundIdx = 0, mNumRounds;
        BetaCircuit mCir;
        span<oc::BetaGate> mGates;

        oc::PRNG mPrng, mPhiPrng;

        span<block> mA, mB, mC, mC2, mD;

        u64 mDebugPrintIdx = -1;
        BetaCircuit::PrintIter mPrint;

        void init(
            u64 n,
            BetaCircuit& cir,
            u64 numThreads,
            u64 batchSize,
            u64 pIdx,
            block seed);

        void setTriples(span<block> a, span<block> b, span<block> c, span<block> d)
        {
            mA = a;
            mB = b;
            mC = c;
            mC2 = c;
            mD = d;
        }

        Proto generateTriple(coproto::Socket& chl);

        template<typename T>
        void setInput(u64 i, oc::MatrixView<T> input)
        {
            static_assert(std::is_trivially_copyable<T>::value, "expecting trivial");
            oc::MatrixView<u8> ii((u8*)input.data(), input.rows(), input.cols() * sizeof(T));
            implSetInput(i, ii, sizeof(T));
        }

        void setZeroInput(u64 i);

        Proto run(coproto::Socket& chl);

        template<typename T>
        void getOutput(u64 i, oc::MatrixView<T> out)
        {
            static_assert(std::is_trivially_copyable<T>::value, "expecting trivial");
            oc::MatrixView<u8> ii((u8*)out.data(), out.rows(), out.cols() * sizeof(T));
            implGetOutput(i, ii, sizeof(T));
        }

        void implSetInput(u64 i, oc::MatrixView<u8> input, u64 alignment);
        void implGetOutput(u64 i, oc::MatrixView<u8> out, u64 alignment);

        oc::MatrixView<u8> getInputView(u64 i);
        oc::MatrixView<u8> getOutputView(u64 i);
        oc::MatrixView<u8> getMemView(BetaBundle& wires);

        SilentTripleGen mSilent;

        u64 numRounds()
        {
            return mNumRounds;
        }

        Proto roundFunction(coproto::Socket& chl);

        Proto multSendP1(span<block> x, coproto::Socket& chl, oc::GateType gt);
        Proto multSendP2(span<block> x, coproto::Socket& chl, oc::GateType gt);

        Proto multRecvP1(span<block> x, span<block> z, coproto::Socket& chl, oc::GateType gt);
        Proto multRecvP2(span<block> x,  span<block> z, coproto::Socket& chl);

        Proto multSend(span<block> x, span<block> y, coproto::Socket& chl, oc::GateType gt)
        {
            if (mIdx == 0)
                return multSendP1(x, y, chl, gt);
            else
                return multSendP2(x, y, chl);
        }
        Proto multSendP1(span<block> x, span<block> y, coproto::Socket& chl, oc::GateType gt);
        Proto multSendP2(span<block> x, span<block> y, coproto::Socket& chl);

        Proto multRecv(span<block> x, span<block> y, span<block> z, coproto::Socket& chl, oc::GateType gt)
        {
            if (mIdx == 0)
                return multRecvP1(x, y, z, chl, gt);
            else
                return multRecvP2(x, y, z, chl);
        }

        Proto multRecvP1(span<block> x, span<block> y, span<block> z, coproto::Socket& chl, oc::GateType gt);
        Proto multRecvP2(span<block> x, span<block> y, span<block> z, coproto::Socket& chl);
    };
}
