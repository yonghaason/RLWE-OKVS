#include "okvs.h"

#include "seal/util/common.h"
#include "seal/util/numth.h"
#include "seal/util/uintarithsmallmod.h"
#include <numeric>

using namespace std;
using namespace seal;
using namespace seal::util;
using namespace oc;

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
        
        // std::cout << "first inner product (before spacing)" << std::endl;
        for (uint32_t i = 0; i < n; ++i)
        {
            const uint64_t *__restrict row = bands_flat.data() + static_cast<size_t>(i) * w;
            uint64_t acc = 0;
            for (uint32_t j = 0; j < w; ++j){
                acc = multiply_add_uint_mod(row[j], x[start_pos[i] + j], acc, modulus);
                // if (i == 0 && j <= 10) {
                //     std::cout << row[j] << " X " << x[start_pos[i] + j] << std::endl;
                // }
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
    const size_t b = m / spacing;

    if (result.size() != n)
        result.resize(n);

    const size_t w = bands_flat.size() / n;

    // std::cout << "first inner product" << std::endl;
    for (uint32_t i = 0; i < n; ++i)
    {
        const uint64_t *__restrict row = bands_flat.data() + static_cast<size_t>(i) * w;
        uint64_t acc = 0;
        auto position = start_pos[i];
            
        for (uint32_t j = 0; j < w; ++j) {
            if (position >= m) position = position - m + 1;
            acc = multiply_add_uint_mod(row[j], x[position], acc, modulus);
            // if (i == 0 && j <= 10) {
            //     std::cout << row[j] << " X " << x[position] << std::endl;
            // }
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
    PRNG prng;
    uint64_t tmp;

    // Generate random-band matrix
#pragma GCC unroll 16 
    for (uint32_t i = 0; i < mN; i++) {
        prng.SetSeed(key[i] ^ seed);
        start_pos[i] = prng.get<uint32_t>() % (mM-mW);
        for (uint32_t j = 0; j < mW; j++) {
            tmp = prng.get<uint64_t>();
            bands_flat[i*mW + j] = barrett_reduce_64(tmp, mModulus);
        }
    }
    setTimePoint("OKVS: Generate Band");
};

    bool PrimeFieldOkvs::sgauss_elimination(
        std::vector<uint64_t> &bands_flat,
        std::vector<uint64_t> &value,
        const std::vector<uint32_t> &start_pos,
        std::vector<uint64_t> &solution)
{   
    // Sort
    vector<uint32_t> perm(mN);
    iota(perm.begin(), perm.end(), 0);
    sort(perm.begin(), perm.end(), [&start_pos](uint32_t a, uint32_t b) 
        { return start_pos[a] < start_pos[b]; });

    vector<const uint64_t *> row_ptr(mN);
    vector<uint64_t> v_perm(mN);
    vector<uint32_t> s_perm(mN);
    for (uint32_t i = 0; i < mN; ++i)
    {
        uint32_t src = perm[i];
        row_ptr[i] = bands_flat.data() + src * mW;
        v_perm[i] = value[src];
        s_perm[i] = start_pos[src];
    }

    vector<uint32_t> offsets(mN);

    for (uint32_t i = 0; i < mN; ++i)
    {
        uint64_t *__restrict row_i = const_cast<uint64_t *>(row_ptr[i]);

        uint32_t offset = UINT32_MAX;
        uint64_t pivot_val = 0;
        for (uint32_t j = 0; j < mW; ++j)
        {
            if (row_i[j])
            {
                pivot_val = row_i[j];
                offset = j;
                break;
            }
        }
        if (offset == UINT32_MAX)
            return false;
        offsets[i] = offset;

        uint64_t inv_val;
        try_invert_uint_mod(pivot_val, mModulus, inv_val);

        row_i[offset] = 1;
        for (uint32_t j = offset + 1; j < mW; ++j)
            row_i[j] = multiply_uint_mod(row_i[j], inv_val, mModulus);
        v_perm[i] = multiply_uint_mod(v_perm[i], inv_val, mModulus);

        uint32_t pivot_col = s_perm[i] + offset;
        for (uint32_t r = i + 1; r < mN; ++r)
        {
            if (s_perm[r] > pivot_col)
                break;

            uint64_t *__restrict row_r = const_cast<uint64_t *>(row_ptr[r]);
            uint32_t idx_r = pivot_col - s_perm[r]; 

            uint64_t factor = row_r[idx_r];
            if (!factor)
                continue;

            row_r[idx_r] = 0;
            for (uint32_t k = pivot_col + 1; k < s_perm[i] + mW; ++k)
            {
                uint32_t i_k = k - s_perm[i];
                uint32_t r_k = k - s_perm[r];

                uint64_t tmp = multiply_uint_mod(row_i[i_k], factor, mModulus);
                row_r[r_k] = sub_uint_mod(row_r[r_k], tmp, mModulus);
            }
            uint64_t tmp_b = multiply_uint_mod(factor, v_perm[i], mModulus);
            v_perm[r] = sub_uint_mod(v_perm[r], tmp_b, mModulus);
        }
    }

    solution.assign(mM, 0);
    for (int i = mN - 1; i >= 0; --i)
    {
        const uint64_t *__restrict row_i = row_ptr[i];
        uint32_t base = s_perm[i];
        uint32_t pivot_idx = base + offsets[i];

        uint64_t acc = 0;
        for (uint32_t k = offsets[i] + 1; k < mW; ++k)
            acc = multiply_add_uint_mod(row_i[k], solution[base + k], acc, mModulus);

        solution[pivot_idx] = sub_uint_mod(v_perm[i], acc, mModulus);
    }
    setTimePoint("OKVS: Solve Lin System");
    return true;
};

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