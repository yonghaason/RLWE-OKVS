// Real-field (Gaussian / full-entropy) RB-OKVS encoding-failure probability vs band width.
//
// For full-entropy coefficients (Gaussian over R, or Z_p over a large field) the encoding
// fails iff the band matrix is rank-deficient, which -- generically -- happens iff the
// sparsity pattern admits no perfect row-to-column matching. For equal-width bands this is an
// interval-matching feasibility test, solved exactly by the leftmost-free greedy (the same
// rule the encoder's triangulation uses). We measure Pr[fail] over random patterns and,
// following Appendix F of the paper, fit lambda = -log2(Pr[fail]) linearly in w and
// extrapolate to lambda = 40.
//
// Self-contained C++; emits CSV (n, eps, w, trials, fails, failprob). Fit/plot in python.

#include <cstdint>
#include <cstdio>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

using std::size_t;

// xoshiro256** seeded by splitmix64
struct Rng {
    uint64_t s[4];
    static uint64_t sm(uint64_t& x){ uint64_t z=(x+=0x9E3779B97F4A7C15ULL);
        z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL; z=(z^(z>>27))*0x94D049BB133111EBULL; return z^(z>>31);}
    void seed(uint64_t k){ for(int i=0;i<4;i++) s[i]=sm(k); }
    static inline uint64_t rotl(uint64_t x,int k){ return (x<<k)|(x>>(64-k)); }
    inline uint64_t next(){ uint64_t r=rotl(s[1]*5,7)*9, t=s[1]<<17;
        s[2]^=s[0]; s[3]^=s[1]; s[1]^=s[2]; s[0]^=s[3]; s[2]^=t; s[3]=rotl(s[3],45); return r; }
    inline uint32_t below(uint32_t b){ return (uint32_t)(((unsigned __int128)(uint32_t)next()*b)>>32); }
};

// versioned union-find "leftmost free point >= x" (no O(m) re-init between trials)
static std::vector<int> nf, stamp;
static int curstamp = 0;
static inline int par(int x){ return stamp[x]==curstamp ? nf[x] : x; }
static inline int findfree(int x){
    int r=x; while(par(r)!=r) r=par(r);
    while(par(x)!=r){ int nx=par(x); nf[x]=r; stamp[x]=curstamp; x=nx; }
    return r;
}
static std::vector<int> starts, cnt;

// one trial: true if encoding FAILS (no perfect interval matching).
// Counting-bucket sweep (O(n+m)) -> process starts in left-to-right order, greedy leftmost
// free point via versioned union-find. No per-trial sort, no O(m) re-init of the UF.
bool one_trial(int n, int m, int w, Rng& rng){
    ++curstamp;
    const int span = m - w + 1;
    starts.resize(n);
    for(int i=0;i<n;i++){ int s=(int)rng.below((uint32_t)span); starts[i]=s; cnt[s]++; }
    bool fail=false;
    for(int s=0; s<span; ++s){
        int c=cnt[s];
        while(c-- > 0){
            int q=findfree(s);
            if(q >= s+w){ fail=true; goto done; }
            nf[q]=q+1; stamp[q]=curstamp;               // use point q
        }
        cnt[s]=0;
    }
done:
    for(int i=0;i<n;i++) cnt[starts[i]]=0;              // O(n) clear (covers any unswept on early-fail)
    return fail;
}

int main(int argc,char**argv){
    std::string csv = argc>1?argv[1]:"failprob.csv";
    FILE* f=std::fopen(csv.c_str(),"w");
    fprintf(f,"n,eps,w,trials,fails,failprob\n");

    const double eps = 1.0;
    struct Cfg{ int logn; long trials; };
    Cfg cfgs[] = { {16, 40000}, {18, 12000}, {20, 3500}, {22, 1000} };  // T*n ~ const

    int m_max=0;
    for(auto&c:cfgs){ long n=1L<<c.logn; int m=(int)llround((1+eps)*n); if(m>m_max)m_max=m; }
    nf.assign(m_max+2,0); stamp.assign(m_max+2,0); cnt.assign(m_max+2,0);

    Rng rng; rng.seed(0xA11CE5);
    for(auto&c:cfgs){
        long n=1L<<c.logn; int m=(int)llround((1+eps)*n);
        for(int w=8; w<=18; ++w){
            long fails=0;
            for(long t=0;t<c.trials;t++) fails += one_trial((int)n,m,w,rng)?1:0;
            double fp=(double)fails/c.trials;
            fprintf(f,"%ld,%.2f,%d,%ld,%ld,%.6g\n",n,eps,w,c.trials,fails,fp); fflush(f);
            printf("n=2^%d w=%-2d  fails=%6ld/%-7ld  failprob=%.4g  lambda=%s\n",
                   c.logn,w,fails,c.trials,fp, fp>0?std::to_string(-log2(fp)).c_str():"inf");
        }
        printf("----\n");
    }
    std::fclose(f);
    printf("wrote %s\n",csv.c_str());
    return 0;
}
