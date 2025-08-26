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

using namespace std;
using namespace oc;
using namespace seal;
using namespace rlweOkvs;

void correct_test(const oc::CLP& cmd)
{
    RpmtSender rpmtSender;
    RpmtReceiver rpmtReceiver;

    oc::Timer timerSender;
    oc::Timer timerReceiver;

    timerSender.setTimePoint("start");
    timerReceiver.setTimePoint("start");
    rpmtSender.setTimer(timerSender);
    rpmtReceiver.setTimer(timerReceiver);
    
    int n = 0;
    auto mNumBatch = 0;
    int logp = 0, numSlots = 0;

    rpmtSender.init(n, logp, numSlots);
    rpmtReceiver.init(n, logp, numSlots);


    // X, Y 랜덤하게 채우고
    // (일정부분은 겹치게)

    vector<block> X(n);
    vector<block> Y(n);

    vector<Plaintext> ptxts;
    rpmtSender.preprocess(Y, ptxts);

    vector<Ciphertext> encoded_in_he(mNumBatch);
    rpmtReceiver.encode_and_encrypt(X, encoded_in_he);

    vector<Ciphertext> decoded_in_he;
    rpmtSender.encrypted_decode(encoded_in_he, ptxts, decoded_in_he);
    
    int their_size = 0;

    oc::BitVector results(their_size); // initialize 할 때 알려주자
    rpmtReceiver.decrypt(decoded_in_he, results);

    // 예상하는 결과랑 result랑 맞는지 확인


    if (cmd.isSet("v")) {
        cout << endl;
        timerSender.setTimePoint("correctness check");
        timerReceiver.setTimePoint("correctness check");
        cout << timerSender << endl;
        cout << timerReceiver << endl;
    } 
}