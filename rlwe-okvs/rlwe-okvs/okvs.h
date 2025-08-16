#pragma once
#include "seal/modulus.h"
#include "seal/randomgen.h"
#include "cryptoTools/Common/block.h"
#include "cryptoTools/Crypto/PRNG.h"

namespace rlweOkvs
{
    class PrimeFieldOkvs
    {
        void init(uint64_t n, uint64_t m, uint32_t w, 
            seal::Modulus &modulus, oc::block seed = oc::ZeroBlock) 
        {
            mN = n;
            mM = m;
            mW = w;
            mPrng.SetSeed(seed);
        };
        
        seal::Modulus mModulus;
        uint64_t mN;
        uint64_t mM;
        uint64_t mW; // random-band weight
        oc::PRNG mPrng;

        public:
            bool encode(
                const std::vector<oc::block> &key,
                const std::vector<uint64_t> &value,
                std::vector<uint64_t> &encoded
            );
            bool decode();
    };
}