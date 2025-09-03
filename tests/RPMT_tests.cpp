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
 
    int n = 1 << 17, logp = 60, numSlots = 1 << 13;

    rpmtSender.init(n, logp, numSlots);
    rpmtReceiver.init(n, logp, numSlots);

    // GPT CODE
    PRNG prng;
    prng.SetSeed(oc::toBlock(0xDEADBEEF12345678ull, 0xBADC0FFEE0DDF00Dull));
    vector<block> X(n);
    vector<block> Y(n);

    prng.get(X.data(), static_cast<u64>(n));
    prng.get(Y.data(), static_cast<u64>(n));

    u64 k = 1 + (prng.get<u64>() % static_cast<u64>(n)); //intersection

    // 서로 다른 인덱스 k개 선택(중복 방지)
    std::unordered_set<size_t> selX, selY;
    selX.reserve(k * 2);
    selY.reserve(k * 2);
    while (selX.size() < k) selX.insert(static_cast<size_t>(prng.get<u64>() % n));
    while (selY.size() < k) selY.insert(static_cast<size_t>(prng.get<u64>() % n));

    // 선택한 k쌍에 대해 Y를 X로 덮어써서 교집합을 만들기
    auto itX = selX.begin();
    auto itY = selY.begin();
    for (; itX != selX.end() && itY != selY.end(); ++itX, ++itY)
        Y[*itY] = X[*itX];
    
    //protocol starts
    timer.setTimePoint("test setting");

    vector<Plaintext> ptxts;
    rpmtSender.preprocess(Y, ptxts);
    timer.setTimePoint("preprocess from sender");

    vector<Ciphertext> encoded_in_he;
    rpmtReceiver.encode_and_encrypt(X, encoded_in_he);
    timer.setTimePoint("encode and encrypt from receiver");

    // 여까기진 정상 작동

    vector<Ciphertext> decoded_in_he;
    rpmtSender.encrypted_decode(encoded_in_he, ptxts, decoded_in_he);
    timer.setTimePoint("encrypted decode from sender");

    oc::BitVector results; // initialize
    rpmtReceiver.decrypt(decoded_in_he, results);
    timer.setTimePoint("decrypt from receiver");

    //GPT CODE

    // 예상하는 결과랑 result랑 맞는지 확인
    // results에서 1의 개수(=히트 수) 카운트
    u64 hit_cnt = 0;
    for (u64 i = 0; i < results.size(); ++i)
        hit_cnt += static_cast<u64>(results[i]);

    // 실제 |X ∩ Y| 계산 (block을 16바이트 문자열 키로 변환해서 set 사용)
    auto keyOf = [](const oc::block& b) -> std::string {
        return std::string(reinterpret_cast<const char*>(&b), sizeof(oc::block));
    };

    std::unordered_set<std::string> setX;
    setX.reserve(static_cast<size_t>(n * 2));
    for (const auto& x : X) setX.insert(keyOf(x));

    u64 inter_cnt = 0;
    for (const auto& y : Y)
        if (setX.find(keyOf(y)) != setX.end()) ++inter_cnt;

    // 숫자 일치 여부 확인
    if (hit_cnt != inter_cnt) {
        std::cerr << "[RPMT test] mismatch: results hits=" << hit_cnt
                  << "  |X∩Y|=" << inter_cnt << std::endl;
        throw RTE_LOC; // 요구사항대로 던짐
    }

    //GPT CODE END

    if (cmd.isSet("v")) {
        cout << endl;
        timer.setTimePoint("correctness check");
        cout << timer << endl;
    } 
}