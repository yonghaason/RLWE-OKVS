// Adapted from the GMW implementation in https://github.com/ladnir/volepsi.
#pragma once

#include "Defines.h"
#include "cryptoTools/Circuit/BetaCircuit.h"

namespace volePSI
{
    using BetaCircuit = oc::BetaCircuit;
    using BetaBundle = oc::BetaBundle;

    BetaCircuit isZeroCircuit(oc::u64 bits);

    BetaCircuit sumThresholdCircuit(oc::u64 bits);

    void isZeroCircuit_Test();
}
