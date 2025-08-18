#pragma once
#include "seal/modulus.h"
#include "seal/randomgen.h"
#include "cryptoTools/Common/block.h"
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Common/Timer.h"

namespace rlweOkvs
{

    void mat_vec_mult_band_flat(
        const std::vector<uint64_t> &bands_flat, 
        const std::vector<uint32_t> &s, 
        const std::vector<uint64_t> &x, 
        const seal::Modulus &modulus,
        std::vector<uint64_t> &result);

    class PrimeFieldOkvs : public oc::TimerAdapter
    {
        seal::Modulus mModulus;
        uint32_t mN;
        uint32_t mM;
        uint32_t mW; // random-band weight

        public:
            void init(uint32_t n, uint32_t m, uint32_t w, seal::Modulus modulus) 
            {
                mN = n;
                mM = m;
                mW = w;
                mModulus = modulus;
            };

            void generate_band(
                const std::vector<oc::block> &key,
                std::vector<uint64_t> &bands_flat,
                std::vector<uint32_t> &start_pos,
                oc::block seed);

            bool sgauss_elimination(
                std::vector<uint64_t> &bands_flat,
                std::vector<uint64_t> &value,
                const std::vector<uint32_t> &start_pos,
                std::vector<uint64_t> &solution);

            bool encode(
                const std::vector<oc::block> &key,
                const std::vector<uint64_t> &value,
                std::vector<uint64_t> &encoded,
                oc::block seed = oc::ZeroBlock);

            void decode(
                const std::vector<oc::block> &key,
                const std::vector<uint64_t> &encoded,
                std::vector<uint64_t> &decoded,
                oc::block seed = oc::ZeroBlock);
    };
}