#include "seal/modulus.h"
#include "seal/randomgen.h"
#include "seal/util/common.h"
#include "seal/util/numth.h"
#include "seal/util/uintarithsmallmod.h"
#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <vector>
#include <filesystem> 
#include <cstdint>
#include <mutex>
#include <sstream>

using namespace std;
using namespace seal;
using namespace seal::util;

// A_flat : 크기 n*L  (row i → &A_flat[i*L])
void gen_random_band_flat(
    std::vector<uint64_t> &A_flat, std::vector<uint64_t> &b, std::vector<int> &s, int n, int m, int L,
    const Modulus &modulus)
{
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> s_dis(0, m - L - 1);
    uint64_t start = (modulus.value() == 2) ? 0 : 1;
    std::uniform_int_distribution<uint64_t> mod_dis(start, modulus.value() - 1);

    A_flat.resize(static_cast<size_t>(n) * L);
    b.resize(n);
    s.resize(n);

    for (int i = 0; i < n; ++i)
    {
        s[i] = static_cast<int>(s_dis(gen));
        uint64_t *row = A_flat.data() + static_cast<size_t>(i) * L;
        for (int j = 0; j < L; ++j)
            row[j] = mod_dis(gen);

        b[i] = mod_dis(gen);
    }
}

// 반환: result = A · x   (A: flat band, 행별 시작열 s, 고정폭 L)
std::vector<uint64_t> mat_vec_mult_band_flat(
    const std::vector<uint64_t> &A_flat, const std::vector<int> &s, const std::vector<uint64_t> &x, int n, int L,
    const Modulus &modulus)
{
    std::vector<uint64_t> result(n, 0);

    for (int i = 0; i < n; ++i)
    {
        const uint64_t *__restrict row = A_flat.data() + static_cast<size_t>(i) * L;
        int base = s[i];

        uint64_t acc = 0;
        for (int j = 0; j < L; ++j)
            acc = multiply_add_uint_mod(row[j], x[base + j], acc, modulus);

        result[i] = acc;
    }
    return result;
}

bool sgauss_modified_flat_ptr(
    std::vector<uint64_t> &A_bands_flat, // flat 행렬
    std::vector<uint64_t> &b, std::vector<int> &s, int L, const Modulus &modulus, std::vector<uint64_t> &x, int m,
    bool debug = false)
{
    /* ------------------------------------------------------------------ */
    /* 0) sort permutation                                               */
    /* ------------------------------------------------------------------ */
    int n = b.size();
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(), [&s](int a, int b) { return s[a] < s[b]; });

    /* ------------------------------------------------------------------ */
    /* 1) row-pointer·b·s 재배열 (복사 X, 포인터 배열만)                  */
    /* ------------------------------------------------------------------ */
    std::vector<const uint64_t *> row_ptr(n); // each -> row base
    std::vector<uint64_t> b_perm(n);
    std::vector<int> s_perm(n);
    for (int i = 0; i < n; ++i)
    {
        int src = perm[i];
        row_ptr[i] = A_bands_flat.data() + src * L; // 포인터만
        b_perm[i] = b[src];
        s_perm[i] = s[src];
    }

    /* ------------------------------------------------------------------ */
    std::vector<int> offsets(n);

    for (int i = 0; i < n; ++i)
    {
        uint64_t *__restrict row_i = const_cast<uint64_t *>(row_ptr[i]);

        /* 2) pivot 찾기 */
        int offset = -1;
        uint64_t pivot_val = 0;
        for (int j = 0; j < L; ++j)
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

        /* 3) invert */
        uint64_t inv_val;
        try_invert_uint_mod(pivot_val, modulus, inv_val);

        /* 4) normalize */
        row_i[offset] = 1;
        for (int j = offset + 1; j < L; ++j)
            row_i[j] = multiply_uint_mod(row_i[j], inv_val, modulus);
        b_perm[i] = multiply_uint_mod(b_perm[i], inv_val, modulus);

        /* 5) eliminate */
        int pivot_col = s_perm[i] + offset;
        for (int r = i + 1; r < n; ++r)
        {
            if (s_perm[r] > pivot_col)
                break;

            uint64_t *__restrict row_r = const_cast<uint64_t *>(row_ptr[r]);
            int idx_r = pivot_col - s_perm[r]; // 위치 in row_r

            uint64_t factor = row_r[idx_r];
            if (!factor)
                continue;

            row_r[idx_r] = 0;
            for (int k = pivot_col + 1; k < s_perm[i] + L; ++k)
            {
                int i_k = k - s_perm[i];
                int r_k = k - s_perm[r];

                uint64_t tmp = multiply_uint_mod(row_i[i_k], factor, modulus);
                row_r[r_k] = sub_uint_mod(row_r[r_k], tmp, modulus);
            }
            uint64_t tmp_b = multiply_uint_mod(factor, b_perm[i], modulus);
            b_perm[r] = sub_uint_mod(b_perm[r], tmp_b, modulus);
        }
    }

    /* ------------------------------------------------------------------ */
    /* 6) back-substitution                                              */
    /* ------------------------------------------------------------------ */
    x.assign(m, 0);
    for (int i = n - 1; i >= 0; --i)
    {
        const uint64_t *__restrict row_i = row_ptr[i];
        int base = s_perm[i];
        int pivot_idx = base + offsets[i];

        uint64_t acc = 0;
        for (int k = offsets[i] + 1; k < L; ++k)
            acc = multiply_add_uint_mod(row_i[k], x[base + k], acc, modulus);

        x[pivot_idx] = sub_uint_mod(b_perm[i], acc, modulus);
    }

    /* ------------------------------------------------------------------ */
    if (debug)
    {
        auto result = mat_vec_mult_band_flat(A_bands_flat, s, x, n, L, modulus);
        for (size_t i = 0; i < n; ++i)
            assert(result[i] == b[i]);
    }
    return true;
}

// struct ResultRecord {
//     int n, logp, m, L, fail_count;
//     uint64_t total_trials;
//     double elapsed_ms;
//     string p_str;
// };

bool sgauss_modified_timer_flat_ptr(
    std::vector<uint64_t> &A_bands_flat, // flat 행렬
    std::vector<uint64_t> &b, std::vector<int> &s, int L, const Modulus &modulus, std::vector<uint64_t> &x, int m,
    std::vector<double> &timer, bool debug = false)
{
    using clk = std::chrono::high_resolution_clock;
    timer.assign(6, 0.0);

    /* ------------------------------------------------------------------ */
    /* 0) sort permutation                                               */
    /* ------------------------------------------------------------------ */
    auto t0 = clk::now();
    int n = b.size();
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(), [&s](int a, int b) { return s[a] < s[b]; });
    timer[0] += std::chrono::duration<double, std::milli>(clk::now() - t0).count();

    /* ------------------------------------------------------------------ */
    /* 1) row-pointer·b·s 재배열 (복사 X, 포인터 배열만)                  */
    /* ------------------------------------------------------------------ */
    t0 = clk::now();
    std::vector<const uint64_t *> row_ptr(n); // each -> row base
    std::vector<uint64_t> b_perm(n);
    std::vector<int> s_perm(n);
    for (int i = 0; i < n; ++i)
    {
        int src = perm[i];
        row_ptr[i] = A_bands_flat.data() + src * L; // 포인터만
        b_perm[i] = b[src];
        s_perm[i] = s[src];
    }
    timer[1] += std::chrono::duration<double, std::milli>(clk::now() - t0).count();

    /* ------------------------------------------------------------------ */
    std::vector<int> offsets(n);

    for (int i = 0; i < n; ++i)
    {
        uint64_t *__restrict row_i = const_cast<uint64_t *>(row_ptr[i]);

        /* 2) pivot 찾기 */
        int offset = -1;
        uint64_t pivot_val = 0;
        for (int j = 0; j < L; ++j)
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

        /* 3) invert */
        t0 = clk::now();
        uint64_t inv_val;
        try_invert_uint_mod(pivot_val, modulus, inv_val);
        timer[2] += std::chrono::duration<double, std::milli>(clk::now() - t0).count();

        /* 4) normalize */
        t0 = clk::now();
        row_i[offset] = 1;
        for (int j = offset + 1; j < L; ++j)
            row_i[j] = multiply_uint_mod(row_i[j], inv_val, modulus);
        b_perm[i] = multiply_uint_mod(b_perm[i], inv_val, modulus);
        timer[3] += std::chrono::duration<double, std::milli>(clk::now() - t0).count();

        /* 5) eliminate */
        t0 = clk::now();
        int pivot_col = s_perm[i] + offset;
        for (int r = i + 1; r < n; ++r)
        {
            if (s_perm[r] > pivot_col)
                break;

            uint64_t *__restrict row_r = const_cast<uint64_t *>(row_ptr[r]);
            int idx_r = pivot_col - s_perm[r]; // 위치 in row_r

            uint64_t factor = row_r[idx_r];
            if (!factor)
                continue;

            row_r[idx_r] = 0;
            for (int k = pivot_col + 1; k < s_perm[i] + L; ++k)
            {
                int i_k = k - s_perm[i];
                int r_k = k - s_perm[r];

                uint64_t tmp = multiply_uint_mod(row_i[i_k], factor, modulus);
                row_r[r_k] = sub_uint_mod(row_r[r_k], tmp, modulus);
            }
            uint64_t tmp_b = multiply_uint_mod(factor, b_perm[i], modulus);
            b_perm[r] = sub_uint_mod(b_perm[r], tmp_b, modulus);
        }
        timer[4] += std::chrono::duration<double, std::milli>(clk::now() - t0).count();
    }

    /* ------------------------------------------------------------------ */
    /* 6) back-substitution                                              */
    /* ------------------------------------------------------------------ */
    t0 = clk::now();
    x.assign(m, 0);
    for (int i = n - 1; i >= 0; --i)
    {
        const uint64_t *__restrict row_i = row_ptr[i];
        int base = s_perm[i];
        int pivot_idx = base + offsets[i];

        uint64_t acc = 0;
        for (int k = offsets[i] + 1; k < L; ++k)
            acc = multiply_add_uint_mod(row_i[k], x[base + k], acc, modulus);

        x[pivot_idx] = sub_uint_mod(b_perm[i], acc, modulus);
    }
    timer[5] += std::chrono::duration<double, std::milli>(clk::now() - t0).count();

    /* ------------------------------------------------------------------ */
    if (debug)
    {
        auto result = mat_vec_mult_band_flat(A_bands_flat, s, x, n, L, modulus);
        for (size_t i = 0; i < n; ++i)
            assert(result[i] == b[i]);
    }
    return true;
}

int many_trials(int logn, double m_ratio, int logp, int trials_total)
{
    /* ----- 파일 이름: "<logn>_<m_ratio>_<logp>" ----- */
    std::ostringstream fname_os;
    fname_os << logn << '_' // 예: 10
             << std::fixed << std::setprecision(1) << m_ratio << '_' // 예: 1.1
             << logp << '_' // 예: 60
             << trials_total; // 예: 100000
    const std::string fname = fname_os.str();

    int n = 1 << logn; // n = 2^{logn}
    int m = static_cast<int>(std::ceil(m_ratio * n));

    /* 이미 파일이 있는지 확인 */
    bool file_exists = std::filesystem::exists(fname);

    /* ios::app 모드 → 이어쓰기.  없으면 새로 생성됨 */
    std::ofstream outfile(fname, std::ios::out | std::ios::app);
    if (!outfile)
    {
        std::cerr << "Cannot open output file \"" << fname << "\"\n";
        return 1;
    }

    /* 새 파일이라면 헤더 한 줄 작성 */
    if (!file_exists)
    {
        outfile << "# L\tlog2_failure_rate\ttrials: " << trials_total << '\n';
    }

    /* ----- 공통 modulus / m 계산 ----- */
    Modulus modulus(2);
    if (logp >= 20)
        modulus = Modulus(PlainModulus::Batching(8192, logp));

    /* ----- 병렬 관련 공통 설정 ----- */
    const unsigned th_cnt = std::thread::hardware_concurrency() ? std::thread::hardware_concurrency() : 8;
    std::cout << "Using " << th_cnt << " threads\n";

    /* ----- L 범위 반복 ----- */
    for (int L = 56; L < 74; L += 3)
    {
        std::cout << "\n[Running n=" << n << ", m=" << m << ", L=" << L << "]\n";

        std::atomic<uint64_t> completed{ 0 };
        std::atomic<uint64_t> fails{ 0 };
        std::mutex ms_mtx;
        double global_ms = 0.0;

        std::mutex io_mtx; // 진행률 출력 락

        /* -------- 스레드 함수 -------- */
        auto worker = [&](unsigned /*tid*/, uint64_t my_trials) {
            double local_ms = 0.0;
            uint64_t local_fail = 0;

            for (uint64_t t = 0; t < my_trials; ++t)
            {
                std::vector<uint64_t> A_flat, b, x;
                std::vector<int> s;

                gen_random_band_flat(A_flat, b, s, n, m, L, modulus);

                auto t0 = std::chrono::high_resolution_clock::now();
                if (!sgauss_modified_flat_ptr(A_flat, b, s, L, modulus, x, m, /*debug=*/false))
                    ++local_fail;
                auto t1 = std::chrono::high_resolution_clock::now();
                local_ms += std::chrono::duration<double, std::milli>(t1 - t0).count();

                completed.fetch_add(1, std::memory_order_relaxed);
            }

            fails.fetch_add(local_fail, std::memory_order_relaxed);

            std::lock_guard<std::mutex> g(ms_mtx);
            global_ms += local_ms;
        };

        /* -------- 진행률 모니터 -------- */
        std::thread monitor([&] {
            uint64_t last_pct = 0;
            while (last_pct < 100)
            {
                uint64_t done = completed.load(std::memory_order_relaxed);
                uint64_t pct = done * 100 / trials_total;
                if (pct > last_pct)
                {
                    last_pct = pct;
                    std::lock_guard<std::mutex> lg(io_mtx);
                    std::cout << "\rProgress: " << pct << "% (" << done << "/" << trials_total << ")" << std::flush;
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(50)); // 폴링 간격
            }
        });

        /* -------- 스레드 실행 -------- */
        uint64_t base = trials_total / th_cnt;
        uint64_t extra = trials_total % th_cnt;
        std::vector<std::thread> pool;
        for (unsigned t = 0; t < th_cnt; ++t)
        {
            uint64_t my = base + (t < extra ? 1 : 0);
            pool.emplace_back(worker, t, my);
        }
        for (auto &th : pool)
            th.join();
        monitor.join();

        /* -------- 결과 집계 -------- */
        double failure_rate = static_cast<double>(fails) / trials_total;
        double lg2_failrate = (failure_rate > 0.0) ? std::log2(failure_rate) : -std::numeric_limits<double>::infinity();

        std::cout << "\n -> Failures: " << fails << " / " << trials_total << " | Failure Rate: " << failure_rate
                  << " = 2^" << lg2_failrate << " | Avg Time per Trial: " << global_ms / trials_total << " ms\n";

        /* -------- 파일 기록 -------- */
        outfile << L << '\t' << lg2_failrate << '\n';
    }

    std::cout << "\nResults written to \"" << fname_os.str() << "\"\n";
    return 0;
}

/* ---------------- 실험 설정 ---------------- */
int main(int argc, char *argv[])
{
    /* ---- 인자 파싱 ---- */
    if (argc < 4)
    {
        std::cerr << "Usage: " << argv[0] << " <logn> <logp> <m_ratio> [trials_total]\n"
                  << "  logn          : log2(n), e.g. 10 for n=1024\n"
                  << "  logp          : plaintext modulus bit-size, e.g. 60\n"
                  << "  m_ratio       : m / n (double)\n"
                  << "  L             : band weight\n"
                  << "  trials_total  : number of trials (optional, default 10000)\n";
        return 1;
    }

    int logn = std::stoi(argv[1]);
    int logp = std::stoi(argv[2]);
    double m_ratio = std::stod(argv[3]);
    int L = std::stoi(argv[4]);
    int trials_total = std::stoi(argv[5]);

    int n = 1 << logn; // n = 2^{logn}
    int m = static_cast<int>(std::ceil(m_ratio * n));

    if (L == 0) many_trials(logn, m_ratio, logp, trials_total);
    else {
        Modulus modulus(2);
        if (logp >= 20)
            modulus = Modulus(PlainModulus::Batching(8192, logp));
        std::vector<uint64_t> A_flat, b, x;
        std::vector<int> s;

        auto stt = std::chrono::high_resolution_clock::now();
        gen_random_band_flat(A_flat, b, s, n, m, L, modulus);
        auto end = std::chrono::high_resolution_clock::now();

        auto time = std::chrono::duration<double, std::milli>(end-stt).count();
        cout << "Gen random band: " << time << " ms" << endl;

        stt = std::chrono::high_resolution_clock::now();
        if (!sgauss_modified_flat_ptr(A_flat, b, s, L, modulus, x, m, false)) 
            cout << "fail" << endl;
        end = std::chrono::high_resolution_clock::now();

        time = std::chrono::duration<double, std::milli>(end-stt).count();
        cout << "Sgauss: " << time << " ms" << endl;
        
    }
}
