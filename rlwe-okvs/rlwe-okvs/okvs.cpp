#include "okvs.h"

#include "seal/util/common.h"
#include "seal/util/numth.h"
#include "seal/util/uintarithsmallmod.h"
#define XXH_INLINE_ALL
#include "xxhash.h"
#include <numeric>

using namespace std;
using namespace seal;
using namespace seal::util;
using namespace oc;


static inline uint64_t invert_uint_mod_prime_u64(uint64_t value, uint64_t modulus)
{
    __int128 t = 0;
    __int128 new_t = 1;

    uint64_t r = modulus;
    uint64_t new_r = value;

    while (new_r != 0)
    {
        const uint64_t q = r / new_r;
        const uint64_t tmp_r = r - q * new_r;

        r = new_r;
        new_r = tmp_r;

        const __int128 tmp_t = t - static_cast<__int128>(q) * new_t;
        t = new_t;
        new_t = tmp_t;
    }

    if (t < 0) {
        t += static_cast<__int128>(modulus);
    }

    return static_cast<uint64_t>(t);
}

namespace rlweOkvs
{
    void mat_vec_mult_band_flat(
        const std::vector<uint64_t> &bands_flat, 
        const std::vector<uint32_t> &start_pos, 
        const std::vector<uint64_t> &x, 
        const seal::Modulus &modulus,
        std::vector<uint64_t> &result)
    {
        auto n = start_pos.size();
        auto w = bands_flat.size() / n;
        
        for (uint32_t i = 0; i < n; ++i)
        {
            const uint64_t *__restrict row = bands_flat.data() + static_cast<size_t>(i) * w;
            uint64_t acc = 0;
            for (uint32_t j = 0; j < w; ++j){
                acc = multiply_add_uint_mod(row[j], x[start_pos[i] + j], acc, modulus);
            }
            result[i] = acc;
        }
    };

    void mat_vec_mult_spaced_band_flat(
        const std::vector<uint64_t> &bands_flat,
        const std::vector<uint32_t> &start_pos,
        uint32_t spacing,
        const seal::Modulus &modulus,
        const std::vector<uint64_t> &x,
        std::vector<uint64_t> &result)
{
    const size_t n = start_pos.size();
    const size_t m = x.size();
    if (m % spacing != 0)
        throw std::invalid_argument("row size should be divisible by t");

    if (result.size() != n)
        result.resize(n);

    const size_t w = bands_flat.size() / n;

    for (uint32_t i = 0; i < n; ++i)
    {
        const uint64_t *__restrict row = bands_flat.data() + static_cast<size_t>(i) * w;
        uint64_t acc = 0;
        auto position = start_pos[i];
            
        for (uint32_t j = 0; j < w; ++j) {
            if (position >= m) position = position - m + 1;
            acc = multiply_add_uint_mod(row[j], x[position], acc, modulus);
            position += spacing;
        }

        result[i] = acc;
    }
}

    void PrimeFieldOkvs::generate_band(
        const std::vector<oc::block> &key,
        std::vector<uint64_t> &bands_flat,
        std::vector<uint32_t> &start_pos,
        oc::block seed)
    {        
        u64 seed_u64 = seed.mData[0];
        // Generate random-band matrix
        
#pragma GCC unroll 16 
        for (uint32_t i = 0; i < mN; ++i) {
            const void* p = key[i].data();
            start_pos[i] = XXH32(p, sizeof(oc::block), seed_u64) % (mM - mW + 1);            
            for (uint32_t j = 0; j < mW; j++) {
                bands_flat[i*mW + j] = barrett_reduce_64(
                    XXH64(p, sizeof(oc::block), seed_u64+j), 
                    mModulus);
            }
        }
    };

bool PrimeFieldOkvs::sgauss_elimination(
    std::vector<uint64_t>& bands_flat,
    std::vector<uint64_t>& value,
    const std::vector<uint32_t>& start_pos,
    std::vector<uint64_t>& solution)
{
    // sort
    std::vector<uint32_t> perm(mN);
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(),
              [&start_pos](uint32_t a, uint32_t b) {
                  return start_pos[a] < start_pos[b];
              });

    // setup
    std::vector<uint64_t*> row_ptr(mN);
    std::vector<uint64_t> v_perm(mN);
    std::vector<uint32_t> s_perm(mN);

    for (uint32_t i = 0; i < mN; ++i) {
        const uint32_t src = perm[i];
        row_ptr[i] = bands_flat.data() + static_cast<size_t>(src) * mW;
        v_perm[i] = value[src];
        s_perm[i] = start_pos[src];
    }

    std::vector<uint32_t> offsets(mN, 0);
    std::vector<uint64_t> inv_perm(mN);

    for (uint32_t i = 0; i < mN; ++i) {
        uint64_t* __restrict row_i = row_ptr[i];

        uint32_t off = offsets[i];
        uint64_t pivot_val = 0;
        while (off < mW) {
            if (row_i[off] != 0) {
                pivot_val = row_i[off];
                break;
            }
            ++off;
        }

        if (off == mW) {
            return false;
        }
        offsets[i] = off;

        const uint64_t inv_val = invert_uint_mod_prime_u64(pivot_val, mModulus.value());
        inv_perm[i] = inv_val;

        const uint32_t pivot_col = s_perm[i] + off;

        for (uint32_t r = i + 1; r < mN; ++r) {
            if (s_perm[r] > pivot_col)
                break;

            uint64_t* __restrict row_r = row_ptr[r];
            const uint32_t idx_r = pivot_col - s_perm[r];
            const uint64_t factor = row_r[idx_r];

            if (!factor)
                continue;

            const uint64_t lambda = multiply_uint_mod(factor, inv_val, mModulus);

            row_r[idx_r] = 0;

            if (offsets[r] <= idx_r)
                offsets[r] = idx_r + 1;

            uint32_t i_k = off + 1;
            uint32_t r_k = idx_r + 1;

            for (; i_k < mW; ++i_k, ++r_k) {
                const uint64_t tmp = multiply_uint_mod(row_i[i_k], lambda, mModulus);
                row_r[r_k] = sub_uint_mod(row_r[r_k], tmp, mModulus);
            }

            const uint64_t tmp_b = multiply_uint_mod(v_perm[i], lambda, mModulus);
            v_perm[r] = sub_uint_mod(v_perm[r], tmp_b, mModulus);
        }
    }

    solution.assign(mM, 0);

    for (int i = static_cast<int>(mN) - 1; i >= 0; --i) {
        
        const uint64_t* __restrict row_i = row_ptr[i];
        const uint32_t base = s_perm[i];
        const uint32_t off = offsets[i];
        const uint32_t pivot_idx = base + off;

        __uint128_t acc128 = 0;

        const uint64_t* __restrict a = row_i + off + 1;
        const uint64_t* __restrict b = solution.data() + base + off + 1;
        const uint64_t* __restrict e = row_i + mW;

        // CAUTION: This accumulating could be problematic when mW is large
        for (; a < e; ++a, ++b) {
            acc128 += static_cast<__uint128_t>(*a) * (*b);
        }

        uint64_t limbs[2] = {
            static_cast<uint64_t>(acc128),
            static_cast<uint64_t>(acc128 >> 64)
        };

        const uint64_t acc = barrett_reduce_128(limbs, mModulus);
        const uint64_t rhs = sub_uint_mod(v_perm[i], acc, mModulus);
        solution[pivot_idx] = multiply_uint_mod(rhs, inv_perm[i], mModulus);
    }

    setTimePoint("OKVS: Solve Lin System");
    return true;
}

    bool PrimeFieldOkvs::encode(
        const vector<block> &key,
        const vector<uint64_t> &value,
        vector<uint64_t> &encoded,
        oc::block seed)
{
    vector<uint64_t> bands_flat(mN*mW);
    vector<uint32_t> start_pos(mN);
    generate_band(key, bands_flat, start_pos, seed);

    vector<uint64_t> value_copy(value);

    return sgauss_elimination(bands_flat, value_copy, start_pos, encoded);
}

    void PrimeFieldOkvs::decode(
        const std::vector<oc::block> &key,
        const std::vector<uint64_t> &encoded,
        std::vector<uint64_t> &decoded,
        oc::block seed)
{
    vector<uint64_t> bands_flat(mN*mW);
    vector<uint32_t> start_pos(mN);
    generate_band(key, bands_flat, start_pos, seed);
    mat_vec_mult_band_flat(bands_flat, start_pos, encoded, mModulus, decoded);
}
}
