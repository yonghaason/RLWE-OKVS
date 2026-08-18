// Adapted from the GMW implementation in https://github.com/ladnir/volepsi.
#include "Circuit.h"
#include <string>

#include "cryptoTools/Circuit/BetaLibrary.h"

using namespace oc;

namespace volePSI
{

    BetaCircuit isZeroCircuit(oc::u64 bits)
    {
        BetaCircuit cd;

        BetaBundle a(bits);

        cd.addInputBundle(a);

        auto ts = [](int s) {return std::to_string(s); };
        u64 step = 1;

        for (u64 i = 0; i < bits; ++i)
            cd.addInvert(a[i]);

        while (step < bits)
        {

            cd.addPrint("\n step " + ts(step) + "\n");
            for (u64 i = 0; i + step < bits; i += step * 2)
            {
                cd.addPrint("a[" + ts(i) + "] & a[" + ts(i + step) + "] -> a[" + ts(i) + "]\n");
                cd.addPrint(a.mWires[i]);
                cd.addPrint(" & ");
                cd.addPrint(a.mWires[i + step]);
                cd.addPrint(" -> ");
                cd.addPrint(a.mWires[i]);

                cd.addGate(a.mWires[i], a.mWires[i + step], oc::GateType::And, a.mWires[i]);
            }

            step *= 2;
        }

        a.mWires.resize(1);
        cd.mOutputs.push_back(a);

        cd.levelByAndDepth();

        return cd;
    }

    BetaCircuit sumThresholdCircuit(oc::u64 bits)
    {
        BetaCircuit cd;

        BetaBundle a(bits), b(bits), t(bits), out(1);
        cd.addInputBundle(a);
        cd.addInputBundle(b);
        cd.addInputBundle(t);
        cd.addOutputBundle(out);

        BetaBundle sum(bits);
        cd.addTempWireBundle(sum);

        BetaBundle noTemps(0);
        oc::BetaLibrary::add_build(cd, a, b, sum, noTemps,
                                   oc::BetaLibrary::IntType::Unsigned,
                                   oc::BetaLibrary::Optimized::Depth);

        oc::BetaLibrary::lessThan_build(cd, t, sum, out,
                                        oc::BetaLibrary::IntType::Unsigned,
                                        oc::BetaLibrary::Optimized::Depth);

        return cd;
    }

    void isZeroCircuit_Test()
    {
        u64 n = 128, tt = 100;
        auto cir = isZeroCircuit(n);

        {
            oc::BitVector bv(n), out(1);
            cir.evaluate({ &bv, 1 }, { &out,1 }, false);

            if (out[0] != 1)
                throw RTE_LOC;
        }

        oc::PRNG prng(oc::ZeroBlock);

        for (u64 i = 0; i < tt; ++i)
        {
            oc::BitVector bv(n), out(1);
            bv.randomize(prng);
            if (bv.hammingWeight() == 0)
                continue;

            cir.evaluate({ &bv, 1 }, { &out,1 }, false);

            if (out[0] != 0)
                throw RTE_LOC;
        }

    }
}
