# Asymptotic corner curves: singlet (M~Bin(k,1/4), cost 3^j C(k,j)) vs Majorana (J~Bin(2k+1,1/2) | parity, u=(J-par)/2,
# cost C(2k+1,j)).  q upper bound = cumulative cost; p exact mass outside best window of d+1 consecutive values.
import math, sys
from math import lgamma, log
def lb(n,j): return lgamma(n+1)-lgamma(j+1)-lgamma(n-j+1)
def lse(xs):
    m=max(xs); return m+math.log(sum(math.exp(x-m) for x in xs))
def best_curve(lp, lc):
    K=len(lp)-1; out=[]
    # windows centred at mode: precompute prefix in log domain via sorting is complex; do direct with pruning
    for d in range(0,K):
        q=lse(lc[:d+1]); best=None
        mode=max(range(K+1),key=lambda i: lp[i])
        for s in range(max(0,mode-d-1),min(K-d,mode+1)+1):
            out_=[x for i,x in enumerate(lp) if i<s or i>s+d]
            if not out_: continue
            e=lse(out_)
            if best is None or e<best: best=e
        if best is None: break
        out.append((q,best))
    return out
k=int(sys.argv[1])
# singlet
lpS=[lb(k,m)+m*log(.25)+(k-m)*log(.75) for m in range(k+1)]
lcS=[lb(k,j)+j*log(3) for j in range(k+1)]
# Majorana: u in 0..k, P(u) = C(2k+1,2u+1)/4^k  (par=1 case)
lpM=[lb(2*k+1,2*u+1)-k*log(4) for u in range(k+1)]
lcM=[lb(2*k+1,j) for j in range(k+1)]
cS=best_curve(lpS,lcS); cM=best_curve(lpM,lcM)
L=k*log(2)
for et in [0.05,0.1,0.2,0.3,0.5,0.8]:
    dS=max([-e/L for q,e in cS if q/L<=et] or [0]); dM=max([-e/L for q,e in cM if q/L<=et] or [0])
    print(f'k={k} eta={et}: delta singlet {dS:.5f} ({dS/et**2:.3f} eta^2)  Majorana {dM:.5f} ({dM/et**2:.3f} eta^2)  ratio {dM/max(dS,1e-12):.2f}')
