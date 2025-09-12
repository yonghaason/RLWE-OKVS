#include "RPMT_tests.h"
#include "rlwe-okvs/rpmt.h"

#include "cryptoTools/Common/Defines.h"
#include "cryptoTools/Common/CLP.h"
#include "cryptoTools/Common/Timer.h"
#include "cryptoTools/Crypto/PRNG.h"

#include "seal/seal.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <unordered_set>

using namespace std;
using namespace oc;
using namespace seal;
using namespace rlweOkvs;

void correct_test(const oc::CLP& cmd)
{
    RpmtSender rpmtSender;
    RpmtReceiver rpmtReceiver;

    oc::Timer timer;

    timer.setTimePoint("start");
    rpmtSender.setTimer(timer);
    rpmtReceiver.setTimer(timer);
 
    int n = 1 << 20, logp = 60, numSlots = 1 << 13;

    PRNG prng;
    prng.SetSeed(oc::ZeroBlock);
    vector<block> X(n);
    vector<block> Y(n);

    prng.get(X.data(), static_cast<u64>(n));
    prng.get(Y.data(), static_cast<u64>(n));

    u64 k = 1 + (prng.get<u64>() % static_cast<u64>(n)); //intersection

    std::unordered_set<size_t> selX, selY;
    selX.reserve(k * 2);
    selY.reserve(k * 2);
    while (selX.size() < k) selX.insert(static_cast<size_t>(prng.get<u64>() % n));
    while (selY.size() < k) selY.insert(static_cast<size_t>(prng.get<u64>() % n));

    auto itX = selX.begin();
    auto itY = selY.begin();
    for (; itX != selX.end() && itY != selY.end(); ++itX, ++itY)
        Y[*itY] = X[*itX];
    
    //protocol starts
    timer.setTimePoint("Input Setting");

    rpmtSender.init(n, logp, numSlots);
    rpmtReceiver.init(n, logp, numSlots);
    

    vector<Plaintext> ptxts_diag;
    vector<Plaintext> ptxts_sdiag;
    rpmtSender.preprocess(Y, ptxts_diag, ptxts_sdiag);
    
    vector<Ciphertext> encoded_in_he;
    vector<Ciphertext> encoded_in_he_shift;
    rpmtReceiver.encode_and_encrypt(X, encoded_in_he, encoded_in_he_shift);
    
    vector<Ciphertext> decoded_in_he;
    rpmtSender.encrypted_decode(encoded_in_he, encoded_in_he_shift,
        ptxts_diag, ptxts_sdiag, decoded_in_he);
    
    oc::BitVector results;
    rpmtReceiver.decrypt(decoded_in_he, results);
    
    u64 real = 0;
    std::unordered_set<oc::block> setX(X.begin(), X.end());
    for (const auto& y : Y)
        if (setX.find(y) != setX.end()) ++real;

    if (results.hammingWeight() != real) {
        cerr << "Protocol =" << results.hammingWeight()
             << " / Expected=" << real << endl;
        throw RTE_LOC; 
    }

    if (cmd.isSet("v")) {
        cout << endl;
        timer.setTimePoint("correctness check");
        cout << timer << endl;
    } 
}