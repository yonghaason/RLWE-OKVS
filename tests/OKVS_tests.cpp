#include "OKVS_Tests.h"
#include "rlwe-okvs/okvs.h"
#include "seal/util/common.h"
#include "seal/util/numth.h"
#include "seal/util/uintarithsmallmod.h"

using namespace oc;
using namespace std;
using namespace seal;
using namespace seal::util;
using namespace rlweOkvs;

void encode_test(const oc::CLP& cmd)
{
    u32 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 10));
    u32 w = cmd.getOr("w", 134);
    u32 m = ceil(cmd.getOr("ratio", 1.16)*n);
    u64 logp = cmd.getOr("logp", 60);
    Modulus p(PlainModulus::Batching(8192, logp));

    oc::Timer timer;
    timer.setTimePoint("start");

    PRNG prng(oc::ZeroBlock);
    vector<block> key(n);
    vector<uint64_t> value(n);

    prng.get<block>(key);
    prng.get<uint64_t>(value);

#pragma GCC unroll 16 
    for (uint32_t i = 0; i < n; i++) {
        value[i] = barrett_reduce_64(value[i], p);
    }

    PrimeFieldOkvs okvs;
    okvs.setTimer(timer);
    okvs.init(n, m, w, p);

    vector<uint64_t> encoded(m);
    vector<uint64_t> bands_flat(n*w);
    vector<uint32_t> start_pos(n);
    okvs.generate_band(key, bands_flat, start_pos, oc::ZeroBlock);

    vector<uint64_t> bands_flat_copy(bands_flat);
    vector<uint64_t> value_copy(value);

    timer.setTimePoint("copy");
    
    okvs.sgauss_elimination(bands_flat_copy, value_copy, start_pos, encoded);

    vector<uint64_t> check(n);
    mat_vec_mult_band_flat(bands_flat, start_pos, encoded, p, check);

    for (int i = 0; i < check.size(); i++) {
        if (check[i] != value[i]) {
            throw RTE_LOC;
        }
    }

    if (cmd.isSet("v")) {
        cout << endl;
        timer.setTimePoint("correctness check");
        cout << timer << endl;
    } 
}

void decode_test(const oc::CLP& cmd)
{
    u32 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 10));
    u32 w = cmd.getOr("w", 134);
    u32 m = ceil(cmd.getOr("ratio", 1.16)*n);
    u64 logp = cmd.getOr("logp", 60);
    Modulus p(PlainModulus::Batching(8192, logp));

    oc::Timer timer;
    timer.setTimePoint("start");

    PRNG prng(oc::ZeroBlock);
    vector<block> key1(n);
    vector<block> key2(n);
    vector<uint64_t> value(n);

    prng.get<block>(key1);
    prng.get<block>(key2);
    prng.get<uint64_t>(value);

    for (uint32_t i = 0; i < n/2; i++) {
        key2[i] = key1[i];
    }

    for (uint32_t i = 0; i < n; i++) {
        value[i] = barrett_reduce_64(value[i], p);
    }

    PrimeFieldOkvs okvs;
    okvs.setTimer(timer);
    okvs.init(n, m, w, p);

    vector<uint64_t> encoded(m);
    okvs.encode(key1, value, encoded);

    vector<uint64_t> decoded(n);
    okvs.decode(key2, encoded, decoded);

    auto cnt = 0;
    for (int i = 0; i < decoded.size(); i++) {
        if (value[i] != decoded[i]) cnt++;
    }
    if (cnt != n/2) throw RTE_LOC;

    if (cmd.isSet("v")) {
        cout << endl;
        cout << timer << endl;
    } 
}
