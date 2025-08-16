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
     
    bool PrimeFieldOkvs::encode(
        const vector<block> &key,
        const vector<uint64_t> &value,
        vector<uint64_t> &encoded
    )
{
    // Generate random-band matrix
    vector<uint64_t> A_bands_flat;    
    vector<uint32_t> starting_idx(mN);

    // Sort
    vector<uint32_t> perm(mN);
    iota(perm.begin(), perm.end(), 0);
    sort(perm.begin(), perm.end(), [&starting_idx](uint32_t a, uint32_t b) 
        { return starting_idx[a] < starting_idx[b]; });

    vector<const uint64_t *> row_ptr(mN); // each -> row base
    vector<uint64_t> v_perm(mN);
    vector<uint32_t> s_perm(mN);
    for (uint32_t i = 0; i < mN; ++i)
    {
        uint32_t src = perm[i];
        row_ptr[i] = A_bands_flat.data() + src * mW; // 포인터만
        v_perm[i] = value[src];
        s_perm[i] = starting_idx[src];
    }

    /* ------------------------------------------------------------------ */
    vector<uint32_t> offsets(mN);

    for (uint32_t i = 0; i < mN; ++i)
    {
        uint64_t *__restrict row_i = const_cast<uint64_t *>(row_ptr[i]);

        uint32_t offset = -1;
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
        if (offset == -1)
            return false;
        offsets[i] = offset;

        uint64_t inv_val;
        try_invert_uint_mod(pivot_val, mModulus, inv_val);

        row_i[offset] = 1;
        for (uint32_t j = offset + 1; j < mW; ++j)
            row_i[j] = multiply_uint_mod(row_i[j], inv_val, mModulus);
        v_perm[i] = multiply_uint_mod(v_perm[i], inv_val, mModulus);

        /* 5) eliminate */
        uint32_t pivot_col = s_perm[i] + offset;
        for (uint32_t r = i + 1; r < mN; ++r)
        {
            if (s_perm[r] > pivot_col)
                break;

            uint64_t *__restrict row_r = const_cast<uint64_t *>(row_ptr[r]);
            uint32_t idx_r = pivot_col - s_perm[r]; // 위치 in row_r

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

    /* ------------------------------------------------------------------ */
    /* 6) back-substitution                                              */
    /* ------------------------------------------------------------------ */
    encoded.assign(mM, 0);
    for (uint32_t i = mN - 1; i >= 0; --i)
    {
        const uint64_t *__restrict row_i = row_ptr[i];
        uint32_t base = s_perm[i];
        uint32_t pivot_idx = base + offsets[i];

        uint64_t acc = 0;
        for (uint32_t k = offsets[i] + 1; k < mW; ++k)
            acc = multiply_add_uint_mod(row_i[k], encoded[base + k], acc, mModulus);

        encoded[pivot_idx] = sub_uint_mod(v_perm[i], acc, mModulus);
    }

    // if (debug)
    // {
    //     auto result = mat_vec_mult_band_flat(A_bands_flat, s, x, mN, mW, mModulus);
    //     for (size_t i = 0; i < mN; ++i)
    //         assert(result[i] == b[i]);
    // }
    return true;
}
}