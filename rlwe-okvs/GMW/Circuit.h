#pragma once
// © 2022 Visa.
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include "Defines.h"
#include "cryptoTools/Circuit/BetaCircuit.h"

namespace volePSI
{
    using BetaCircuit = oc::BetaCircuit;
    using BetaBundle = oc::BetaBundle;

    BetaCircuit isZeroCircuit(oc::u64 bits);

    // Threshold comparison over arithmetic shares. Three input bundles of
    // `bits` wires each -- a, b and t -- and a single output wire set to
    // 1(a + b > t), where a + b is computed modulo 2^bits. The two parties
    // feed their arithmetic shares of the cardinality into a and b (each one
    // sets its own bundle and zeroes the other, so the XOR-shared bundle
    // values are the shares themselves) and the public threshold into t.
    // Testing "at least T" therefore means passing t = T - 1 (T >= 1); the
    // strict form is what keeps the output wire free of an invert flag, which
    // the GMW backend would drop.
    BetaCircuit sumThresholdCircuit(oc::u64 bits);

    void isZeroCircuit_Test();
}