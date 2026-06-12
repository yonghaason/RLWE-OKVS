// CKKS multi-block N-spaced (RSB) homomorphic decode, rotation-key variant.
//
// Extends ckks_decode.cpp to m = b*N (b blocks): p is packed into b CKKS ciphertexts and the
// RSB column permutation pi_N spaces each band by N. A band element of the query with
// contiguous start s touches  block (s+j) mod b , wrap-level r = (s%b + j) div b , at output
// slot tau = s div b ; crossing a block boundary shifts the slot by +1 (the super-diagonal),
// which is exactly the homomorphic rotation rot_r. Decode of a layer:
//     out = sum_{blk, r} PlainMult( rot_r(c_blk), diag[blk][r] )
// the rot_r(c_blk) being computed ONCE by Galois-key rotation (not pre-rotated by R) and
// reused across all layers. Queries (= R's keys; S holds the same elements, so every query
// matches) are sequenced into layers by tau so each layer has <= 1 query per slot.
//
// Verifies v'_i == v_i for ALL queries. Garbage removal still deferred.

#include <cstdio>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include "seal/seal.h"
using std::vector; using namespace seal;

// ---- contiguous Gaussian band + min-norm encode (banded Cholesky), starts may collide ----
struct Band { int m,n,w; vector<int> start; vector<double> coeff; };
Band gen_band(int m,int n,int w,std::mt19937_64&rng){
    Band B; B.m=m;B.n=n;B.w=w; B.start.resize(n); B.coeff.assign((size_t)n*w,0.0);
    std::uniform_int_distribution<int> sd(0,m-w); std::normal_distribution<double> g(0,1);
    vector<int> s(n); for(int i=0;i<n;i++) s[i]=sd(rng); std::sort(s.begin(),s.end());
    B.start=s; for(int i=0;i<n;i++)for(int j=0;j<w;j++) B.coeff[(size_t)i*w+j]=g(rng);
    return B;
}
bool minnorm(const Band&B,const vector<double>&v,vector<double>&p){
    int n=B.n,w=B.w; int kd=0;
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){ if(B.start[j]>=B.start[i]+w)break; kd=std::max(kd,j-i);}
    vector<double> A((size_t)n*(kd+1),0.0);
    auto At=[&](int i,int d)->double&{return A[(size_t)i*(kd+1)+d];};
    for(int i=0;i<n;i++)for(int d=0;d<=kd&&d<=i;d++){int j=i-d;
        int lo=std::max(B.start[i],B.start[j]),hi=std::min(B.start[i]+w,B.start[j]+w); double a=0;
        for(int c=lo;c<hi;c++) a+=B.coeff[(size_t)i*w+(c-B.start[i])]*B.coeff[(size_t)j*w+(c-B.start[j])];
        At(i,d)=a; }
    for(int j=0;j<n;j++){ double s=At(j,0); for(int d=1;d<=kd&&d<=j;d++){double x=At(j,d);s-=x*x;}
        if(s<=0)return false; double L=std::sqrt(s); At(j,0)=L;
        for(int i=j+1;i<=j+kd&&i<n;i++){int dij=i-j; double t=At(i,dij);
            for(int d=1;d<=kd;d++){int k=j-d;if(k<0)break;int di=i-k,dj=j-k;if(di<=kd&&dj<=kd)t-=At(i,di)*At(j,dj);}
            At(i,dij)=t/L; } }
    vector<double> y(v);
    for(int i=0;i<n;i++){double s=y[i];for(int d=1;d<=kd&&d<=i;d++)s-=At(i,d)*y[i-d];y[i]=s/At(i,0);}
    for(int i=n-1;i>=0;i--){double s=y[i];for(int d=1;d<=kd&&i+d<n;d++)s-=At(i+d,d)*y[i+d];y[i]=s/At(i,0);}
    p.assign(B.m,0.0); for(int i=0;i<n;i++){int st=B.start[i];const double*r=&B.coeff[(size_t)i*w];
        for(int j=0;j<w;j++) p[st+j]+=r[j]*y[i]; }
    return true;
}

int main(){
    const int LOGRING=13; const size_t poly=1<<LOGRING; const int N=poly/2; // N=4096 slots
    const int b=8;                  // blocks
    const int m=b*N;                // = 32768
    const int n=m/2;                // eps=1  -> n=16384
    const int w=32;                 // band width (>> 2^-40 threshold at this n)
    const int maxr=(b-1 + w-1)/b;   // max wrap level
    std::mt19937_64 rng(2024);

    Band B=gen_band(m,n,w,rng);
    std::normal_distribution<double> g(0,1);
    vector<double> v(n); for(auto&x:v)x=g(rng);
    vector<double> p; if(!minnorm(B,v,p)){printf("encode failed (rank); raise w\n");return 1;}
    double nv=0,np=0; for(double x:v)nv+=x*x; for(double x:p)np+=x*x; nv=sqrt(nv);np=sqrt(np);
    printf("encode: n=%d m=%d=b*N (b=%d N=%d) w=%d eps=%.2f  ||p||/||v||=%.4f  maxr=%d\n",
           n,m,b,N,w,(double)(m-n)/n, np/nv, maxr);

    // pi_N permutation: block blk holds  phat_blk[slot] = p[slot*b + blk]
    vector<vector<double>> phat(b, vector<double>(N,0.0));
    for(int blk=0;blk<b;blk++) for(int slot=0;slot<N;slot++) phat[blk][slot]=p[slot*b+blk];

    // sequencing by tau = start div b  -> layers (<=1 query per slot per layer)
    // diag[layer][blk][r] : length-N plaintext diagonal
    vector<int> slotuse(N,0); int L=0;
    struct Cell{ int layer,slot; };
    vector<Cell> qcell(n);
    for(int i=0;i<n;i++){ int tau=B.start[i]/b; qcell[i]={slotuse[tau]++, tau}; L=std::max(L,qcell[i].layer+1); }
    printf("sequencing: %d layers (n/N=%.1f)\n", L, (double)n/N);
    // diag[L][b][maxr+1] each size N
    vector<vector<vector<vector<double>>>> diag(L, vector<vector<vector<double>>>(b,
        vector<vector<double>>(maxr+1, vector<double>(N,0.0))));
    for(int i=0;i<n;i++){ int s=B.start[i], tau=s/b, ell=qcell[i].layer; const double* c=&B.coeff[(size_t)i*w];
        for(int j=0;j<w;j++){ int blk=(s+j)%b, r=(s%b + j)/b; diag[ell][blk][r][tau]+=c[j]; } }

    // ---- CKKS ----
    EncryptionParameters parms(scheme_type::ckks); parms.set_poly_modulus_degree(poly);
    parms.set_coeff_modulus(CoeffModulus::Create(poly,{60,40,40,60}));
    SEALContext context(parms); KeyGenerator keygen(context);
    SecretKey sk=keygen.secret_key(); PublicKey pk; keygen.create_public_key(pk);
    GaloisKeys gk; keygen.create_galois_keys(gk);
    Encryptor enc(context,pk); Evaluator ev(context); Decryptor dec(context,sk); CKKSEncoder ecd(context);
    double scale=pow(2.0,40); size_t slots=ecd.slot_count();
    printf("CKKS: ring=%zu slots=%zu scale=2^40\n",poly,slots);

    // base blocks
    vector<Ciphertext> c(b);
    for(int blk=0;blk<b;blk++){ vector<double> blkv(slots,0.0); for(int t=0;t<N;t++) blkv[t]=phat[blk][t];
        Plaintext pt; ecd.encode(blkv,scale,pt); enc.encrypt(pt,c[blk]); }
    // homomorphic rotations rot_r(c_blk), r=1..maxr  (computed once, reused across layers)
    vector<vector<Ciphertext>> rc(b, vector<Ciphertext>(maxr+1));
    for(int blk=0;blk<b;blk++){ rc[blk][0]=c[blk];
        for(int r=1;r<=maxr;r++) ev.rotate_vector(c[blk], r, gk, rc[blk][r]); }

    // decode each layer:  out = sum_{blk,r} PlainMult(rot_r(c_blk), diag[ell][blk][r])
    vector<vector<double>> res(L);
    for(int ell=0;ell<L;ell++){ Ciphertext acc; bool first=true;
        for(int blk=0;blk<b;blk++) for(int r=0;r<=maxr;r++){
            const vector<double>& d=diag[ell][blk][r];
            bool nz=false; for(double x:d){ if(x!=0.0){nz=true;break;} } if(!nz) continue;
            vector<double> dd(slots,0.0); for(int t=0;t<N;t++) dd[t]=d[t];
            Plaintext pt; ecd.encode(dd, rc[blk][r].scale(), pt); ev.mod_switch_to_inplace(pt, rc[blk][r].parms_id());
            Ciphertext term; ev.multiply_plain(rc[blk][r], pt, term); ev.rescale_to_next_inplace(term);
            term.scale()=scale;
            if(first){acc=term;first=false;} else ev.add_inplace(acc,term);
        }
        Plaintext ptr; dec.decrypt(acc,ptr); ecd.decode(ptr,res[ell]);
    }

    // verify v'_i = res[layer_i][slot_i] == v_i for all queries
    double maxe=0,rmse=0; for(int i=0;i<n;i++){ double e=fabs(res[qcell[i].layer][qcell[i].slot]-v[i]);
        maxe=std::max(maxe,e); rmse+=e*e; } rmse=sqrt(rmse/n);
    printf("decode (multi-block, rotation-key): max|v'-v|=%.3e  rms|v'-v|=%.3e  (~%.1f bits)\n",
           maxe,rmse,-log2(rmse/(nv/sqrt(n))));
    printf("  (%d base ctxts + %d homomorphic rotations, %d layers)\n", b, b*maxr, L);
    return 0;
}
