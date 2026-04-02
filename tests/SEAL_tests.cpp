#include "SEAL_tests.h"
#include "rlwe-okvs/okvs.h"
#include "seal/seal.h"
#include "seal/util/common.h"
#include "seal/util/numth.h"
#include "seal/util/uintarithsmallmod.h"

using namespace oc;
using namespace std;
using namespace seal;
using namespace seal::util;
using namespace rlweOkvs;

void context_consistency_test(const oc::CLP& cmd)
{
    u32 nslots = cmd.getOr("n", 1ull << cmd.getOr("nn", 12));    
    u64 logp = cmd.getOr("logp", 40);

    vector<int> coeffs{50, 30, 20, 20};

    EncryptionParameters parms(scheme_type::bgv);
    parms.set_poly_modulus_degree(nslots);
    parms.set_coeff_modulus(CoeffModulus::Create(nslots, coeffs));
    parms.set_plain_modulus(PlainModulus::Batching(nslots, logp));
    SEALContext context(parms, true, seal::sec_level_type::none);

    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();

    BatchEncoder encoder(context);
    Encryptor encryptor(context, secret_key);
    Decryptor decryptor(context, secret_key);

    Plaintext ptxt;
    vector<uint64_t> plain_vec(nslots);
    iota(plain_vec.begin(), plain_vec.end(), 0);
    encoder.encode(plain_vec, ptxt);
    Ciphertext ctxt;
    encryptor.encrypt_symmetric(ptxt, ctxt);

    // Below are assumed to be performed by the other party
    EncryptionParameters parms2(scheme_type::bgv);
    parms2.set_poly_modulus_degree(nslots);
    parms2.set_coeff_modulus(CoeffModulus::Create(nslots, coeffs));
    parms2.set_plain_modulus(PlainModulus::Batching(nslots, logp));
    SEALContext context2(parms2, true, seal::sec_level_type::none);

    BatchEncoder encoder2(context2);
    Evaluator evaluator(context2);
    
    vector<uint64_t> plain_vec2(nslots);
    Plaintext ptxt2;
    iota(plain_vec2.begin(), plain_vec2.end(), 1);    
    encoder2.encode(plain_vec2, ptxt2);
    Ciphertext ctxt_mult; 
    evaluator.multiply_plain(ctxt, ptxt2, ctxt_mult);

    Ciphertext ctxt_mult3(ctxt_mult);
    evaluator.mod_switch_to_next_inplace(ctxt_mult3);
    evaluator.mod_switch_to_next_inplace(ctxt_mult3);
   
    // Now the secret-key holder decrypt the result.    
    Plaintext dec;
    vector<uint64_t> dec_vec;
    decryptor.decrypt(ctxt_mult, dec);
    encoder.decode(dec, dec_vec);
    for (size_t i = 0; i < dec_vec.size(); i++) {
        if (dec_vec[i] != plain_vec[i] * plain_vec2[i])
            throw RTE_LOC;
    }
        
    decryptor.decrypt(ctxt_mult3, dec);
    encoder.decode(dec, dec_vec);
    for (size_t i = 0; i < dec_vec.size(); i++) {
        if (dec_vec[i] != plain_vec[i] * plain_vec2[i])
            throw RTE_LOC;
    }
    if (cmd.isSet("v")) {
        cout << endl;
        cout << "Noise budgets" << endl;
        cout << " - Fresh: " 
            << decryptor.invariant_noise_budget(ctxt) << endl;
        cout << " - " << ctxt_mult.coeff_modulus_size() << "-th lv: " 
            << decryptor.invariant_noise_budget(ctxt_mult) << endl;
        cout << " - " << ctxt_mult3.coeff_modulus_size() << "-th lv: " 
            << decryptor.invariant_noise_budget(ctxt_mult3) << endl;
    }
};