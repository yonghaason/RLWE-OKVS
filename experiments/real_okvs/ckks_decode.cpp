// CKKS homomorphic decode of a real-field (Gaussian) RB-OKVS, rotation-key variant.
//
// Demonstrates the encode -> CKKS encrypt (scale 2^40) -> homomorphic spaced/diagonalized
// band decode -> decrypt round-trip that delivers  <row(y), p> = payload  in the encrypted
// state, using HOMOMORPHIC ROTATIONS (Galois keys) -- NOT pre-rotated ciphertexts.
//
// Scope of this first cut (kept faithful to the RSB idea, simplifications noted):
//   * single block: m = N = #CKKS slots, p in one ciphertext;
//   * single diagonalized layer: the n queries occupy distinct diagonal positions
//     pos_i in [0, N-w], so the decode is  c* = sum_{j=0}^{w-1} rot_j(c_p) (*) diag_j ,
//     i.e. w-1 homomorphic rotations + w plain-mults (one multiplicative level);
//   * Gaussian width-w bands, min-norm encoding (banded Cholesky).
// Deferred to next steps (see design note): collisions -> sequencing into layers; m>N ->
// N-spacing across blocks; and garbage removal (indicator OKVS + ctxt x ctxt mask).
//
// Build (see run command in REPORT); links SEAL 4.1 + HEXL.

#include <cstdint>
#include <cstdio>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include "seal/seal.h"

using std::vector;
using namespace seal;

// ---------- min-norm real-field band OKVS encode (banded Cholesky) ----------
// Row i: support [pos_i, pos_i+w), Gaussian coeffs; pos_i distinct in [0, N-w] (diagonalized).
// Solve min-norm  p* = M^T (MM^T)^{-1} v  in R^N.
struct Band { int N, n, w; vector<int> pos; vector<double> coeff; }; // coeff[i*w+j]

Band gen_band(int N, int n, int w, std::mt19937_64& rng){
    Band B; B.N=N; B.n=n; B.w=w; B.pos.resize(n); B.coeff.assign((size_t)n*w,0.0);
    // distinct positions in [0, N-w]
    vector<int> all(N-w+1); std::iota(all.begin(),all.end(),0);
    std::shuffle(all.begin(),all.end(),rng);
    for(int i=0;i<n;i++) B.pos[i]=all[i];
    std::sort(B.pos.begin(),B.pos.end());
    std::normal_distribution<double> g(0.0,1.0);
    for(int i=0;i<n;i++) for(int j=0;j<w;j++) B.coeff[(size_t)i*w+j]=g(rng);
    return B;
}

// banded symmetric (lower) storage
static bool chol_solve_minnorm(const Band& B, const vector<double>& v, vector<double>& p){
    int n=B.n,w=B.w;
    // bandwidth of MM^T (rows sorted by pos): max index gap with |pos_i-pos_j|<w
    int kd=0; for(int i=0;i<n;i++){ for(int j=i+1;j<n;j++){ if(B.pos[j]>=B.pos[i]+w) break; kd=std::max(kd,j-i);} }
    // build A=MM^T (lower banded): A[i][i-d]
    vector<double> A((size_t)n*(kd+1),0.0);
    auto At=[&](int i,int d)->double&{ return A[(size_t)i*(kd+1)+d]; };
    for(int i=0;i<n;i++){ for(int d=0;d<=kd&&d<=i;d++){ int j=i-d;
        int lo=std::max(B.pos[i],B.pos[j]), hi=std::min(B.pos[i]+w,B.pos[j]+w);
        double acc=0; for(int c=lo;c<hi;c++) acc+=B.coeff[(size_t)i*w+(c-B.pos[i])]*B.coeff[(size_t)j*w+(c-B.pos[j])];
        At(i,d)=acc; } }
    // banded Cholesky A=LL^T (in place)
    for(int j=0;j<n;j++){ double s=At(j,0); for(int d=1;d<=kd&&d<=j;d++){double x=At(j,d); s-=x*x;}
        if(s<=0) return false; double Ljj=std::sqrt(s); At(j,0)=Ljj;
        for(int i=j+1;i<=j+kd&&i<n;i++){ int dij=i-j; double t=At(i,dij);
            for(int d=1;d<=kd;d++){ int k=j-d; if(k<0)break; int di=i-k,dj=j-k; if(di<=kd&&dj<=kd) t-=At(i,di)*At(j,dj);}
            At(i,dij)=t/Ljj; } }
    // solve LL^T y = v
    vector<double> y(v);
    for(int i=0;i<n;i++){ double s=y[i]; for(int d=1;d<=kd&&d<=i;d++) s-=At(i,d)*y[i-d]; y[i]=s/At(i,0);}
    for(int i=n-1;i>=0;i--){ double s=y[i]; for(int d=1;d<=kd&&i+d<n;d++) s-=At(i+d,d)*y[i+d]; y[i]=s/At(i,0);}
    // p = M^T y
    p.assign(B.N,0.0);
    for(int i=0;i<n;i++){ int s=B.pos[i]; const double* row=&B.coeff[(size_t)i*w]; for(int j=0;j<w;j++) p[s+j]+=row[j]*y[i]; }
    return true;
}

int main(){
    // ---- instance ----
    const int LOGRING=14; const size_t poly=1<<LOGRING;     // ring dim 2^14 -> N=8192 slots
    const int N=poly/2;                                       // CKKS slots = encoding length m
    const int w=24;
    const int n=4000;                                         // queries; eps=(N-n)/n approx 1.05
    std::mt19937_64 rng(12345);

    Band B=gen_band(N,n,w,rng);
    std::normal_distribution<double> g(0.0,1.0);
    vector<double> v(n); for(auto&x:v) x=g(rng);              // real payloads ~ N(0,1)

    vector<double> p;
    if(!chol_solve_minnorm(B,v,p)){ printf("encode failed (rank)\n"); return 1; }
    double nv=0,np=0,pinf=0; for(double x:v)nv+=x*x; for(double x:p){np+=x*x;pinf=std::max(pinf,std::fabs(x));}
    nv=std::sqrt(nv); np=std::sqrt(np);
    printf("encode: n=%d N(=m=slots)=%d w=%d eps=%.3f  ||p||/||v||=%.4f  ||p||_inf=%.4f\n",
           n,N,w,(double)(N-n)/n, np/nv, pinf);

    // ---- CKKS setup, scale 2^40 ----
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly);
    parms.set_coeff_modulus(CoeffModulus::Create(poly, {60,40,40,60}));
    SEALContext context(parms);
    KeyGenerator keygen(context);
    SecretKey sk=keygen.secret_key();
    PublicKey pk; keygen.create_public_key(pk);
    GaloisKeys gk; keygen.create_galois_keys(gk);             // rotation (Galois) keys
    Encryptor encryptor(context,pk); Evaluator evaluator(context); Decryptor decryptor(context,sk);
    CKKSEncoder encoder(context);
    double scale=std::pow(2.0,40);
    size_t slots=encoder.slot_count();
    printf("CKKS: ring=%zu slots=%zu scale=2^40 coeffmod_bits={60,40,40,60}\n",poly,slots);

    // encrypt p (one block, N slots)
    vector<double> pblk(slots,0.0); for(int i=0;i<N;i++) pblk[i]=p[i];
    Plaintext ptp; encoder.encode(pblk,scale,ptp);
    Ciphertext cp; encryptor.encrypt(ptp,cp);

    // ---- rotation-key decode: c* = sum_{j=0}^{w-1} rot_j(cp) (*) diag_j ----
    // diag_j[pos_i] = coeff_{i,j}, else 0.
    Ciphertext acc; bool first=true;
    for(int j=0;j<w;j++){
        Ciphertext rot;
        if(j==0) rot=cp; else evaluator.rotate_vector(cp,j,gk,rot);   // homomorphic rotation
        vector<double> diag(slots,0.0);
        for(int i=0;i<n;i++) diag[B.pos[i]]=B.coeff[(size_t)i*w+j];
        Plaintext ptd; encoder.encode(diag, rot.scale(), ptd);
        evaluator.mod_switch_to_inplace(ptd, rot.parms_id());
        Ciphertext term; evaluator.multiply_plain(rot,ptd,term);
        evaluator.rescale_to_next_inplace(term);
        term.scale()=scale;                                          // force exact scale match
        if(first){ acc=term; first=false; } else evaluator.add_inplace(acc,term);
    }

    // ---- decrypt and compare v'_i = result[pos_i] vs v_i ----
    Plaintext ptr; decryptor.decrypt(acc,ptr);
    vector<double> res; encoder.decode(ptr,res);
    double maxerr=0,rmserr=0; for(int i=0;i<n;i++){ double e=std::fabs(res[B.pos[i]]-v[i]); maxerr=std::max(maxerr,e); rmserr+=e*e; }
    rmserr=std::sqrt(rmserr/n);
    printf("decode (rotation-key): max|v'-v|=%.3e  rms|v'-v|=%.3e  (=> ~%.1f bits preserved)\n",
           maxerr,rmserr, -std::log2(rmserr/ (nv/std::sqrt(n)) ));
    printf("predicted Corollary 1 bound  sqrt(w)*eta with eta~||p||_inf*2^-40 : %.3e\n",
           std::sqrt((double)w)*pinf*std::pow(2.0,-40));
    return 0;
}
