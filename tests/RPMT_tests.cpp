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

namespace {
    // rpmt.cpp과 같은 방식으로 m과 배치 수를 계산
    inline uint32_t compute_m(uint32_t n, uint64_t numSlots) {
        // m = ((ceil(1.16 * n) / numSlots) + 1) * numSlots;
        return static_cast<uint32_t>(((std::ceil(1.16 * n) / numSlots) + 1) * numSlots);
    }
    inline uint32_t compute_num_batch(uint32_t m, uint64_t numSlots) {
        return static_cast<uint32_t>(m / numSlots);
    }

    // 고정 seed 로 블록 벡터 채우기
    inline void fill_random_blocks(std::vector<block>& v, const block& seed = oc::ZeroBlock) {
        PRNG prng(seed);
        prng.get<block>(v.data(), v.size());
    }
}

// -----------------------------------------------------------------------------
// Receiver-side smoke test
//   - RpmtReceiver::init()
//   - RpmtReceiver::encode_and_encrypt()
//   산출된 ciphertext 개수가 mNumBatch와 동일한지 검증
// -----------------------------------------------------------------------------
void receiver_test(const oc::CLP& cmd)
{
    const uint32_t n     = cmd.getOr("n", 1ull << cmd.getOr("nn", 10)); // 기본 2^10
    const uint64_t logp  = cmd.getOr("logp", 60ull);
    const uint64_t slots = cmd.getOr("slots", 8192ull);

    // 입력 집합 X 준비
    std::vector<block> X(n);
    fill_random_blocks(X, oc::ZeroBlock);

    // 수신자 초기화 및 인코딩+암호화
    RpmtReceiver recv;
    recv.init(n, static_cast<uint32_t>(logp), slots);

    std::vector<Ciphertext> ctxts;
    recv.encode_and_encrypt(X, ctxts);

    // 기대 개수 확인
    const uint32_t m         = compute_m(n, slots);
    const uint32_t num_batch = compute_num_batch(m, slots);

    if (ctxts.size() != num_batch) {
        std::cout << "[RPMT Receiver Test] Expected " << num_batch
                  << " ciphertexts but got " << ctxts.size() << std::endl;
        throw RTE_LOC;
    }

    if (cmd.isSet("v")) {
        std::cout << "[RPMT Receiver Test]\n"
                  << " n=" << n << " logp=" << logp << " slots=" << slots << "\n"
                  << " m=" << m << " batches=" << num_batch << "\n"
                  << " -> encode_and_encrypt() produced " << ctxts.size() << " ciphertexts\n";
    }
}

// -----------------------------------------------------------------------------
// Sender-side integration test (복호화 없이 파이프라인 확인)
//   파이프라인: Receiver::encode_and_encrypt -> Sender::preprocess -> Sender::encrypted_decode
//   체크: decoded_in_he.size() == L (= ptxts.size() / num_batch)
// -----------------------------------------------------------------------------
void sender_test(const oc::CLP& cmd)
{
    const uint32_t n     = cmd.getOr("n", 1ull << cmd.getOr("nn", 10));
    const uint64_t logp  = cmd.getOr("logp", 60ull);
    const uint64_t slots = cmd.getOr("slots", 8192ull);

    // X, Y를 부분 교집합으로 생성 (앞 절반은 동일하게)
    std::vector<block> X(n), Y(n);
    {
        PRNG prng1(oc::ZeroBlock);
        PRNG prng2(oc::toBlock(123456789ull, 987654321ull));
        prng1.get<block>(X.data(), X.size());
        prng2.get<block>(Y.data(), Y.size());
        for (uint32_t i = 0; i < n / 2; ++i) {
            Y[i] = X[i];
        }
    }

    // 양측 초기화
    RpmtReceiver recv;
    RpmtSender   send;
    recv.init(n, static_cast<uint32_t>(logp), slots);
    send.init(n, static_cast<uint32_t>(logp), slots);

    // Receiver: 인코딩 + 암호문 생성
    std::vector<Ciphertext> encoded_in_he;
    recv.encode_and_encrypt(X, encoded_in_he);

    // Sender: 미리 연산할 plaintext 배치 생성
    std::vector<Plaintext> ptxts;
    send.preprocess(Y, ptxts);

    // Sender: 암호문 × 평문곱 + 덧셈으로 디코딩된 암호문 벡터 생성
    std::vector<Ciphertext> decoded_in_he;
    send.encrypted_decode(encoded_in_he, ptxts, decoded_in_he);

    // 개수 검증: L = ptxts.size() / num_batch
    const uint32_t m         = compute_m(n, slots);
    const uint32_t num_batch = compute_num_batch(m, slots);
    const uint32_t L         = static_cast<uint32_t>(ptxts.size() / num_batch);

    if (decoded_in_he.size() != L) {
        std::cout << "[RPMT Sender Test] Expected " << L
                  << " decoded ciphertexts but got " << decoded_in_he.size() << std::endl;
        throw RTE_LOC;
    }

    if (cmd.isSet("v")) {
        std::cout << "[RPMT Sender Test]\n"
                  << " n=" << n << " logp=" << logp << " slots=" << slots << "\n"
                  << " m=" << m << " batches=" << num_batch << " L=" << L << "\n"
                  << " -> preprocess() produced " << ptxts.size() << " plaintexts\n"
                  << " -> encrypted_decode() produced " << decoded_in_he.size() << " ciphertexts\n";
    }
}
