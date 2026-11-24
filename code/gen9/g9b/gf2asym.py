"""GF(2) one-bit scheme: row rank of A_i = rank of e_{t-i}(N_1..N_t) on (J2+2J1)^{(x)t}
 = sum_n C(t,n) 2^{t-n} r(n,d), r(n,d)=rank of mult by e_d in F2[N_1..N_n]/(N_s^2) = sum_s rank_2 W_{s,s+d}(n).
Wilson 1990: rank_2 W_{s,k}(n) = sum_{j<=s, C(k-j,s-j) odd} C(n,j)-C(n,j-1)   (s<=min(k,n-k))."""
from math import comb, log2
import itertools, sys
from functools import lru_cache
def gf2rank(rows):
    r=0; rows=list(rows)
    basis={}
    for v in rows:
        while v:
            h=v.bit_length()-1
            if h in basis: v^=basis[h]
            else: basis[h]=v; r+=1; break
    return r
def W_brute(n,s,k):
    S=[sum(1<<i for i in c) for c in itertools.combinations(range(n),s)]
    K=[sum(1<<i for i in c) for c in itertools.combinations(range(n),k)]
    rows=[]
    for a in S:
        v=0
        for j,b in enumerate(K):
            if a&b==a: v|=1<<j
        rows.append(v)
    return gf2rank(rows)
@lru_cache(None)
def W_wilson(n,s,k):
    if k>n or s>k: return 0
    if s>n-k: return W_wilson(n,n-k,n-s)
    return sum(comb(n,j)-(comb(n,j-1) if j>0 else 0) for j in range(s+1) if comb(k-j,s-j)%2==1)
def WQ(n,s,k):
    if k>n or s>k: return 0
    return min(comb(n,s),comb(n,k))
if __name__=="__main__":
    bad=0
    for n in range(1,11):
        for s in range(n+1):
            for k in range(s,n+1):
                if W_brute(n,s,k)!=W_wilson(n,s,k): bad+=1
    print("Wilson formula vs brute force, n<=10: mismatches",bad)
    def r(n,d): return sum(W_wilson(n,s,s+d) for s in range(0,n-d+1))
    def rowrank2(t,i):
        d=t-i; return sum(comb(t,n)*2**(t-n)*r(n,d) for n in range(t+1))
    def rowrankQ(t,i): return sum(comb(t,z)*3**z for z in range(i+1))
    print("t=5 GF2 row ranks:",[rowrank2(5,i) for i in range(6)], "(direct: 1,12,66,232,496,1024)")
    print("  t   rho_Q(A_i)  rho_GF2(A_i)  argmax_i/t")
    for t in (2,3,4,5,8,16,32,64,128,256,512,1024):
        cq=max(log2(rowrankQ(t,i))+log2(rowrankQ(t,t-i)) for i in range(t+1))/(4*t)
        vals=[log2(rowrank2(t,i))+log2(rowrank2(t,t-i)) for i in range(t+1)]
        c2=max(vals)/(4*t); print("%5d  %.4f     %.4f       %.3f"%(t,cq,c2,vals.index(max(vals))/t)); sys.stdout.flush()
