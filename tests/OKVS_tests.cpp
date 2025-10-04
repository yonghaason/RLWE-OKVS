#include "OKVS_Tests.h"
#include "rlwe-okvs/okvs.h"
#include "seal/util/common.h"
#include "seal/util/numth.h"
#include "seal/util/uintarithsmallmod.h"


#include <mutex>
#include <filesystem> 
#include <fstream>
#include <random>

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

    for (u64 i = 0; i < check.size(); i++) {
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

    u32 cnt = 0;
    for (u64 i = 0; i < decoded.size(); i++) {
        if (value[i] != decoded[i]) cnt++;
    }
    if (cnt != n/2) throw RTE_LOC;

    if (cmd.isSet("v")) {
        cout << endl;
        cout << timer << endl;
    } 
}

void width_test(const oc::CLP& cmd)
{
    u32 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 10));
    double m_ratio = cmd.getOr("ratio", 1.16);
    u64 logp = cmd.getOr("logp", 60);
    u32 trials_total = cmd.getOr("trials", 100000);
    u32 w0 = cmd.getOr("w0", 0);
    u32 w1 = cmd.getOr("w1", 0);

    if (w0 == w1) w1 = w0 + 1;
    if (w1 == 0) w1 = w0 + 1;

    std::ostringstream fname_os;
    fname_os << log2ceil(n) << '_' 
             << std::fixed << std::setprecision(1) << m_ratio << '_' 
             << logp << '_'
             << trials_total; 
    const std::string fname = fname_os.str();

    u32 m = ceil(m_ratio*n/8192)*8192;  // t_s = 2^13 

    cout << "m : " << m << endl;

    bool file_exists = std::filesystem::exists(fname);

    std::ofstream outfile(fname, std::ios::out | std::ios::app);
    if (!outfile)
    {
        std::cerr << "Cannot open output file \"" << fname << "\"\n";
        throw RTE_LOC;
    }

    if (!file_exists)
    {
        outfile << "# L\tlog2_failure_rate\ttrials: " << trials_total << '\n';
    }

    Modulus p(2);
    if (logp != 1)
        p = Modulus(PlainModulus::Batching(1024, logp));

    const unsigned th_cnt = std::thread::hardware_concurrency() ? std::thread::hardware_concurrency() : 8;
    std::cout << "Using " << th_cnt << " threads\n";

    oc::PRNG prng0(oc::ZeroBlock);

    for (u32 w = w0; w < w1; w++)
    {
        std::cout << "\n[Running n=" << n << ", m=" << m << ", w=" << w << "]\n";

        std::atomic<uint64_t> completed{ 0 };
        std::atomic<uint64_t> fails{ 0 };
        std::mutex ms_mtx;
        double global_ms = 0.0;

        std::mutex io_mtx;

        auto worker = [&](unsigned /*tid*/, uint64_t my_trials) {
            double local_ms = 0.0;
            uint64_t local_fail = 0;

            oc::PRNG prng(prng0.get());

            for (uint64_t t = 0; t < my_trials; ++t)
            {
                vector<u64> A_flat(n*w), b(n);
                vector<u32> s(n);

                PrimeFieldOkvs okvs;
                okvs.init(n, m, w, p);

                prng.get<u64>(A_flat);
                prng.get<u64>(b);
                prng.get<u32>(s);

                for (uint32_t i = 0; i < A_flat.size(); i++) {
                    A_flat[i] = barrett_reduce_64(A_flat[i], p);
                }

                for (uint32_t i = 0; i < n; i++) {
                    b[i] = barrett_reduce_64(b[i], p);
                    s[i] %= m-w-1;
                }

                vector<u64> x;
                
                auto t0 = std::chrono::high_resolution_clock::now();
                if (!okvs.sgauss_elimination(A_flat, b, s, x))
                    ++local_fail;

                auto t1 = std::chrono::high_resolution_clock::now();
                local_ms += std::chrono::duration<double, std::milli>(t1 - t0).count();

                completed.fetch_add(1, std::memory_order_relaxed);
            }

            fails.fetch_add(local_fail, std::memory_order_relaxed);

            std::lock_guard<std::mutex> g(ms_mtx);
            global_ms += local_ms;
        };

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
                std::this_thread::sleep_for(std::chrono::milliseconds(50)); 
            }
        });

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

        double failure_rate = static_cast<double>(fails) / trials_total;
        double lg2_failrate = (failure_rate > 0.0) ? std::log2(failure_rate) : -std::numeric_limits<double>::infinity();

        std::cout << "\n -> Failures: " << fails << " / " << trials_total << " | Failure Rate: " << failure_rate
                  << " = 2^" << lg2_failrate << " | Avg Time per Trial: " << global_ms / trials_total << " ms\n";

        outfile << w << '\t' << lg2_failrate << '\n';
    }

    std::cout << "\nResults written to \"" << fname_os.str() << "\"\n";
}