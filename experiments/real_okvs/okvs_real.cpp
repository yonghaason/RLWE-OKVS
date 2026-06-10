// Real-field RB-OKVS conditioning experiment.
//
// This is a faithful port of rlwe-okvs/rlwe-okvs/okvs.cpp (PrimeFieldOkvs) from the
// prime field Z_p to the real field R (double), so that "how large does the encoding
// vector get?" becomes a meaningful question. We compare:
//   (1) the naive triangular solve  -- the exact algorithm of sgauss_elimination,
//       ported to double (free variables set to 0), and
//   (2) the minimum-l2-norm solve   -- p* = M^T (M M^T)^{-1} v via banded Cholesky,
//       i.e. the O(n w^2) method enabled by M M^T being banded after sorting by start.
// We also measure sigma_min(M), sigma_max(M), and the garbage |<row(y'),p>| on random
// non-key queries (the quantity CKKS masking must keep small).
//
// Self-contained: depends only on the C++ standard library. No SEAL / Eigen.
//
// Band conventions mirror PrimeFieldOkvs::generate_band:
//   - start in [0, m-w], width w, NO wraparound (plain RB matrix; the RSB column
//     permutation preserves singular values, so it does not change any norm here).
//   - BINARY band: w bits, leading bit forced to 1 (anchors the start), rest ~ U{0,1}.
//   - GAUSS  band: w real coeffs ~ N(0,1) (real analogue of full-field Z_p coeffs).

#include <cstdint>
#include <cstdio>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <string>

using std::vector;
using std::size_t;

enum class BandType { Binary, Gauss };

struct BandMatrix {
    int n;            // rows (#keys)
    int m;            // cols (encoding length)
    int w;            // band width
    vector<int> start;        // start[i]      in [0, m-w]
    vector<double> coeff;     // coeff[i*w + j] band coefficient at column start[i]+j
};

// ---- band generation (port of generate_band; bandtype-configurable) ----
BandMatrix gen_band(int n, int m, int w, BandType bt, std::mt19937_64& rng) {
    BandMatrix M; M.n = n; M.m = m; M.w = w;
    M.start.resize(n);
    M.coeff.assign((size_t)n * w, 0.0);
    std::uniform_int_distribution<int> startD(0, m - w);
    std::uniform_int_distribution<int> bit(0, 1);
    std::normal_distribution<double> gauss(0.0, 1.0);
    for (int i = 0; i < n; ++i) {
        M.start[i] = startD(rng);
        for (int j = 0; j < w; ++j) {
            double c;
            if (bt == BandType::Binary) c = (j == 0) ? 1.0 : (double)bit(rng);
            else                        c = gauss(rng);
            M.coeff[(size_t)i * w + j] = c;
        }
    }
    return M;
}

// ---- decode: result[i] = <row_i, x> ----
double dot_row(const BandMatrix& M, int i, const vector<double>& x) {
    const double* row = &M.coeff[(size_t)i * M.w];
    int s = M.start[i];
    double acc = 0.0;
    for (int j = 0; j < M.w; ++j) acc += row[j] * x[s + j];
    return acc;
}

// ---- naive triangular solve over double: faithful port of sgauss_elimination ----
// Returns false if a row has no available pivot (encoding failure). Free cols stay 0.
bool naive_triangular(const BandMatrix& Min, const vector<double>& value,
                      vector<double>& solution) {
    const int n = Min.n, m = Min.m, w = Min.w;
    vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(),
              [&](int a, int b){ return Min.start[a] < Min.start[b]; });

    // working copies (sorted)
    vector<double> rows((size_t)n * w);
    vector<double> v_perm(n);
    vector<int>    s_perm(n);
    for (int i = 0; i < n; ++i) {
        int src = perm[i];
        std::copy(&Min.coeff[(size_t)src*w], &Min.coeff[(size_t)src*w]+w, &rows[(size_t)i*w]);
        v_perm[i] = value[src];
        s_perm[i] = Min.start[src];
    }
    vector<int> offsets(n, 0);
    vector<double> inv_pivot(n, 0.0);

    const double EPS = 1e-300; // exact-zero test (pivots are generically nonzero over R)
    for (int i = 0; i < n; ++i) {
        double* row_i = &rows[(size_t)i*w];
        int off = offsets[i];
        double pivot_val = 0.0;
        while (off < w) { if (std::fabs(row_i[off]) > EPS) { pivot_val = row_i[off]; break; } ++off; }
        if (off == w) return false;
        offsets[i] = off;
        inv_pivot[i] = 1.0 / pivot_val;
        int pivot_col = s_perm[i] + off;

        for (int r = i + 1; r < n; ++r) {
            if (s_perm[r] > pivot_col) break;            // outside band window
            double* row_r = &rows[(size_t)r*w];
            int idx_r = pivot_col - s_perm[r];
            if (idx_r < 0 || idx_r >= w) continue;
            double factor = row_r[idx_r];
            if (factor == 0.0) continue;
            double lambda = factor * inv_pivot[i];
            row_r[idx_r] = 0.0;
            if (offsets[r] <= idx_r) offsets[r] = idx_r + 1;
            int ik = off + 1, rk = idx_r + 1;
            for (; ik < w; ++ik, ++rk) row_r[rk] -= row_i[ik] * lambda;
            v_perm[r] -= v_perm[i] * lambda;
        }
    }

    solution.assign(m, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        const double* row_i = &rows[(size_t)i*w];
        int base = s_perm[i], off = offsets[i];
        int pivot_idx = base + off;
        double acc = 0.0;
        for (int k = off + 1; k < w; ++k) acc += row_i[k] * solution[base + k];
        solution[pivot_idx] = (v_perm[i] - acc) * inv_pivot[i];
    }
    return true;
}

// ---- symmetric banded matrix (lower), LAPACK-like: A[i][d] = full(i, i-d), d=0..kd ----
struct SymBanded {
    int n, kd;
    vector<double> a; // size n*(kd+1); a[i*(kd+1)+d]
    double& at(int i, int d) { return a[(size_t)i*(kd+1)+d]; }
    double  at(int i, int d) const { return a[(size_t)i*(kd+1)+d]; }
};

// Build (perm-sorted) M M^T in banded form. Returns the row permutation used (sorted by
// start) so we can map v and recover p in original column order.
SymBanded build_MMt(const BandMatrix& M, const vector<int>& perm, int kd) {
    int n = M.n;
    SymBanded A; A.n = n; A.kd = kd; A.a.assign((size_t)n*(kd+1), 0.0);
    // For each pair (i>=j) in sorted order with band overlap, dot product of the two rows.
    for (int i = 0; i < n; ++i) {
        int pi = perm[i], si = M.start[pi];
        const double* ri = &M.coeff[(size_t)pi*M.w];
        for (int d = 0; d <= kd && d <= i; ++d) {
            int j = i - d, pj = perm[j], sj = M.start[pj];
            // overlap of [si,si+w) and [sj,sj+w):
            int lo = std::max(si, sj), hi = std::min(si + M.w, sj + M.w);
            if (lo >= hi) { A.at(i,d) = 0.0; continue; }
            const double* rj = &M.coeff[(size_t)pj*M.w];
            double acc = 0.0;
            for (int c = lo; c < hi; ++c) acc += ri[c - si] * rj[c - sj];
            A.at(i,d) = acc;
        }
    }
    return A;
}

// In-place banded Cholesky: A = L L^T, L stored in same lower-banded layout. false if not SPD.
bool chol_banded(SymBanded& A) {
    int n = A.n, kd = A.kd;
    for (int j = 0; j < n; ++j) {
        double sum = A.at(j,0);
        for (int d = 1; d <= kd && d <= j; ++d) { double v = A.at(j,d); sum -= v*v; }
        if (sum <= 0.0) return false;
        double Ljj = std::sqrt(sum);
        A.at(j,0) = Ljj;
        for (int i = j+1; i <= j+kd && i < n; ++i) {
            int dij = i - j;
            double s = A.at(i,dij);
            for (int d = 1; d <= kd; ++d) {
                int k = j - d; if (k < 0) break;
                int di = i - k, dj = j - k;
                if (di <= kd && dj <= kd) s -= A.at(i,di) * A.at(j,dj);
            }
            A.at(i,dij) = s / Ljj;
        }
    }
    return true;
}

// Solve A x = b given Cholesky factor L in A (lower banded). b overwritten / x returned.
void chol_solve(const SymBanded& L, vector<double>& b) {
    int n = L.n, kd = L.kd;
    // forward: L y = b
    for (int i = 0; i < n; ++i) {
        double s = b[i];
        for (int d = 1; d <= kd && d <= i; ++d) s -= L.at(i,d) * b[i-d];
        b[i] = s / L.at(i,0);
    }
    // backward: L^T x = y
    for (int i = n-1; i >= 0; --i) {
        double s = b[i];
        for (int d = 1; d <= kd && i+d < n; ++d) s -= L.at(i+d,d) * b[i+d];
        b[i] = s / L.at(i,0);
    }
}

// Symmetric-banded mat-vec  y = A x  (A given by ORIGINAL banded entries, not factor).
void symband_matvec(const SymBanded& A, const vector<double>& x, vector<double>& y) {
    int n = A.n, kd = A.kd;
    std::fill(y.begin(), y.end(), 0.0);
    for (int i = 0; i < n; ++i) {
        double s = A.at(i,0) * x[i];
        for (int d = 1; d <= kd && d <= i; ++d) {
            s += A.at(i,d) * x[i-d];      // lower
            y[i-d] += A.at(i,d) * x[i];   // upper (symmetric)
        }
        y[i] += s;
    }
}

static double norm2(const vector<double>& x){ double s=0; for(double v:x) s+=v*v; return std::sqrt(s); }
static double norminf(const vector<double>& x){ double s=0; for(double v:x) s=std::max(s,std::fabs(v)); return s; }

struct Metrics {
    double ok;                 // naive triangulation success (1/0)
    double ratio_naive, ratio_min;     // ||p||2 / ||v||2
    double linf_naive, linf_min;       // ||p||inf
    double sig_min, sig_max, cond;
    double bound_min;          // ||v||/sigma_min   (theory upper bound on ||p*||)
    double res_naive, res_min; // max|M p - v|
    double g_naive_rms, g_naive_max, g_min_rms, g_min_max; // garbage on non-keys
};

Metrics run_once(int n, int m, int w, BandType bt, int n_garbage, std::mt19937_64& rng) {
    Metrics R{};
    BandMatrix M = gen_band(n, m, w, bt, rng);
    std::normal_distribution<double> gauss(0.0,1.0);
    vector<double> v(n); for (auto& x : v) x = gauss(rng);
    double nv = norm2(v);

    // (1) naive triangular
    vector<double> p_naive;
    bool ok = naive_triangular(M, v, p_naive);
    R.ok = ok ? 1.0 : 0.0;
    if (!ok) return R;

    // sorted-by-start permutation and MMt bandwidth (max index gap among overlapping rows)
    vector<int> perm(n); std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(), [&](int a,int b){ return M.start[a] < M.start[b]; });
    int kd = 0;
    for (int i = 0; i < n; ++i)
        for (int j = i+1; j < n; ++j) {
            if (M.start[perm[j]] >= M.start[perm[i]] + w) break;
            kd = std::max(kd, j - i);
        }

    // (2) min-norm: solve (MMt) y = v_perm, p* = M^T y
    SymBanded A = build_MMt(M, perm, kd);
    SymBanded L = A;                 // copy for factorization
    bool spd = chol_banded(L);
    vector<double> p_min(m, 0.0);
    if (spd) {
        vector<double> y(n);
        for (int i = 0; i < n; ++i) y[i] = v[perm[i]];
        chol_solve(L, y);            // y = (MMt)^{-1} v_perm
        for (int i = 0; i < n; ++i) { // p* = M^T y
            int pi = perm[i], s = M.start[pi];
            const double* row = &M.coeff[(size_t)pi*w];
            for (int j = 0; j < w; ++j) p_min[s+j] += row[j] * y[i];
        }
    }

    // sigma_max(M)^2 = lambda_max(MMt)  via power iteration
    // sigma_min(M)^2 = lambda_min(MMt)  via inverse power iteration (uses L)
    auto rayleigh = [&](const vector<double>& x){
        vector<double> Ax(n); symband_matvec(A, x, Ax);
        double num=0,den=0; for(int i=0;i<n;++i){num+=x[i]*Ax[i];den+=x[i]*x[i];}
        return num/den;
    };
    vector<double> z(n); for (auto& x : z) x = gauss(rng);
    double lmax = 0;
    for (int it = 0; it < 300; ++it) {
        vector<double> Az(n); symband_matvec(A, z, Az);
        double nz = norm2(Az); if (nz==0) break; for(int i=0;i<n;++i) z[i]=Az[i]/nz;
        lmax = rayleigh(z);
    }
    double lmin = 0;
    if (spd) {
        vector<double> q(n); for (auto& x : q) x = gauss(rng);
        for (int it = 0; it < 300; ++it) {
            vector<double> t = q; chol_solve(L, t);   // t = (MMt)^{-1} q
            double nt = norm2(t); if (nt==0) break; for(int i=0;i<n;++i) q[i]=t[i]/nt;
        }
        lmin = rayleigh(q);
    }
    R.sig_max = std::sqrt(std::max(lmax,0.0));
    R.sig_min = std::sqrt(std::max(lmin,0.0));
    R.cond = (R.sig_min>0)? R.sig_max/R.sig_min : INFINITY;
    R.bound_min = (R.sig_min>0)? nv / R.sig_min : INFINITY;

    R.ratio_naive = norm2(p_naive)/nv;
    R.ratio_min   = norm2(p_min)/nv;
    R.linf_naive  = norminf(p_naive);
    R.linf_min    = norminf(p_min);

    // decode residuals (sanity: should be ~0)
    double rn=0, rm=0;
    for (int i=0;i<n;++i){ rn=std::max(rn,std::fabs(dot_row(M,i,p_naive)-v[i]));
                           rm=std::max(rm,std::fabs(dot_row(M,i,p_min)-v[i])); }
    R.res_naive=rn; R.res_min=rm;

    // garbage on random non-key queries: fresh random bands, evaluate <row', p>
    std::uniform_int_distribution<int> startD(0, m - w);
    std::uniform_int_distribution<int> bit(0,1);
    double gn2=0,gnmax=0,gm2=0,gmmax=0;
    for (int q=0;q<n_garbage;++q){
        int s=startD(rng); double an=0, am=0;
        for(int j=0;j<w;++j){ double c=(bt==BandType::Binary)?((j==0)?1.0:(double)bit(rng)):gauss(rng);
                              an+=c*p_naive[s+j]; am+=c*p_min[s+j]; }
        gn2+=an*an; gnmax=std::max(gnmax,std::fabs(an));
        gm2+=am*am; gmmax=std::max(gmmax,std::fabs(am));
    }
    R.g_naive_rms=std::sqrt(gn2/n_garbage); R.g_naive_max=gnmax;
    R.g_min_rms  =std::sqrt(gm2/n_garbage); R.g_min_max =gmmax;
    return R;
}

// average metrics over reps
Metrics run_avg(int n, double eps, int w, BandType bt, int reps, int n_garbage, uint64_t seed){
    int m = (int)std::llround((1.0+eps)*n);
    Metrics S{}; int good=0;
    std::mt19937_64 rng(seed);
    double ok=0;
    Metrics acc{};
    for (int r=0;r<reps;++r){
        Metrics M = run_once(n,m,w,bt,n_garbage,rng);
        ok += M.ok;
        if (M.ok < 0.5) continue;
        good++;
        acc.ratio_naive+=M.ratio_naive; acc.ratio_min+=M.ratio_min;
        acc.linf_naive +=M.linf_naive;  acc.linf_min +=M.linf_min;
        acc.sig_min+=M.sig_min; acc.sig_max+=M.sig_max; acc.cond+=M.cond;
        acc.bound_min+=M.bound_min; acc.res_naive+=M.res_naive; acc.res_min+=M.res_min;
        acc.g_naive_rms+=M.g_naive_rms; acc.g_naive_max+=M.g_naive_max;
        acc.g_min_rms+=M.g_min_rms;     acc.g_min_max+=M.g_min_max;
    }
    S.ok = ok/reps;
    if (good>0){ double g=good;
        S.ratio_naive=acc.ratio_naive/g; S.ratio_min=acc.ratio_min/g;
        S.linf_naive=acc.linf_naive/g;   S.linf_min=acc.linf_min/g;
        S.sig_min=acc.sig_min/g; S.sig_max=acc.sig_max/g; S.cond=acc.cond/g;
        S.bound_min=acc.bound_min/g; S.res_naive=acc.res_naive/g; S.res_min=acc.res_min/g;
        S.g_naive_rms=acc.g_naive_rms/g; S.g_naive_max=acc.g_naive_max/g;
        S.g_min_rms=acc.g_min_rms/g;     S.g_min_max=acc.g_min_max/g;
    }
    return S;
}

const char* bname(BandType b){ return b==BandType::Binary?"binary":"gauss"; }

int main(int argc, char** argv){
    std::string csv = (argc>1)? argv[1] : "results.csv";
    FILE* f = std::fopen(csv.c_str(), "w");
    fprintf(f, "sweep,band,n,eps,w,m,ok,ratio_naive,ratio_min,linf_naive,linf_min,"
               "sig_min,sig_max,cond,bound_min,res_naive,res_min,"
               "g_naive_rms,g_naive_max,g_min_rms,g_min_max\n");
    auto emit = [&](const char* sweep, BandType bt, int n, double eps, int w, const Metrics& M){
        int m=(int)std::llround((1.0+eps)*n);
        fprintf(f,"%s,%s,%d,%.4f,%d,%d,%.4f,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.3g,%.3g,%.6g,%.6g,%.6g,%.6g\n",
                sweep,bname(bt),n,eps,w,m,M.ok,M.ratio_naive,M.ratio_min,M.linf_naive,M.linf_min,
                M.sig_min,M.sig_max,M.cond,M.bound_min,M.res_naive,M.res_min,
                M.g_naive_rms,M.g_naive_max,M.g_min_rms,M.g_min_max);
        printf("[%-7s %-6s n=%-5d eps=%.2f w=%-3d] ok=%.2f  ||p||/||v||: naive=%9.3g  min=%7.3f  "
               "sig_min=%8.2e  cond=%9.2e  garbage(min)=%.3g\n",
               sweep,bname(bt),n,eps,w,M.ok,M.ratio_naive,M.ratio_min,M.sig_min,M.cond,M.g_min_rms);
        fflush(f);
    };

    const int REPS = 10, NG = 256;
    uint64_t seed = 0xC0FFEE;

    // Sweep A: vary w (binary, n=1024, eps=1.0)
    for (int w : {10,12,14,16,20,24,30,40,50})
        emit("w", BandType::Binary, 1024, 1.0, w, run_avg(1024,1.0,w,BandType::Binary,REPS,NG,seed++));

    // Sweep B: vary eps (binary, n=1024, w=24)
    for (double eps : {0.05,0.1,0.2,0.3,0.5,1.0,2.0})
        emit("eps", BandType::Binary, 1024, eps, 24, run_avg(1024,eps,24,BandType::Binary,REPS,NG,seed++));

    // Sweep C: vary n (binary, eps=1.0, w=24)
    for (int n : {256,512,1024,2048,4096})
        emit("n", BandType::Binary, n, 1.0, 24, run_avg(n,1.0,24,BandType::Binary,REPS,NG,seed++));

    // Sweep D: binary vs gauss bands (n=1024, eps=1.0)
    for (int w : {12,16,24,32}) {
        emit("band", BandType::Binary, 1024, 1.0, w, run_avg(1024,1.0,w,BandType::Binary,REPS,NG,seed++));
        emit("band", BandType::Gauss,  1024, 1.0, w, run_avg(1024,1.0,w,BandType::Gauss, REPS,NG,seed++));
    }

    std::fclose(f);
    printf("\nWrote %s\n", csv.c_str());
    return 0;
}
