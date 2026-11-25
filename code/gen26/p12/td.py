# Tensor domination test: Smith(Per_k(E_{a,i})) vs union_{k0+k1=k} Smith(Per_k0(E_{a-1,i-1})) (+) Smith(Per_k1(E_{a-1,i-1}))
import sys,itertools
from fractions import Fraction as Fr
from lib12 import *
from levelm3 import perZ
NB=300;MOD=1<<NB
def perprofile(E,kmax):
    Ez=[[toZ(z,MOD) for z in r] for r in E]; n=len(E); out={0:[0]}
    for k in range(1,kmax+1):
        subs=list(itertools.combinations(range(n),k))
        out[k]=smith_vals([[perZ([[Ez[r][u] for u in U] for r in R]) for U in subs] for R in subs],NB)
    return out
a=int(sys.argv[1]); t=int(sys.argv[2]); zeta=Fr(sys.argv[3]); kmax=int(sys.argv[4])
E,c,T,B=cauchy_block(a,t,zeta,0); n=len(E)
E0,c0,_,_=cauchy_block(a-1,t//2,zeta,0)
# scalar: E0 entries valuation align: v2(c0) vs v2(c)-1
sh=v2q(c)-1-v2q(c0)
big=perprofile(E,kmax); small=perprofile(E0,kmax)
for k in range(1,kmax+1):
    tens=sorted(x+y+k*sh for k0 in range(k+1) for x in small.get(k0,[]) for y in small.get(k-k0,[]) if k0<=len(E0) and k-k0<=len(E0))
    sp=big[k]
    le=all(p<=q for p,q in zip(sp,tens)); ge=all(p>=q for p,q in zip(sp,tens))
    print(f"a={a} t={t} k={k}: Per_k <= tensor: {le}  >= : {ge}  equal: {sp==tens}")
    if not le: print("   big ",sp[:40],"\n   tens",tens[:40])
