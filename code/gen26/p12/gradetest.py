# Direct test of the grading consequence of Thm U: for random F in Lambda(M) (integer comb. of rows of M), the degree-k part of its
# Moebius transform lies in rowspace(Per_k(E)) (lattice membership via Smith index).
import sys,itertools,random
from fractions import Fraction as Fr
from lib12 import *
from fastA import smith_fast, NB, MOD
random.seed(3)
def mob(vec,n):
    v=vec[:]
    for j in range(n):
        for X in range(1<<n):
            if X>>j&1: v[X]-=v[X^(1<<j)]
    return v
def member(rows,v):
    s1=smith_fast(rows); s2=smith_fast(rows+[v])[:len(s1)]
    return s1==s2
def test(E,label,trials=3):
    n=len(E); M=Mmat(E); Mz=[[toZ(z,MOD) for z in r] for r in M]
    ok=True
    for tr in range(trials):
        coef=[random.randint(-3,3) for _ in range(1<<n)]
        F=[sum(cf*Mz[S][X] for S,cf in enumerate(coef))%MOD for X in range(1<<n)]
        Fm=mob(F,n)
        for k in range(1,n+1):
            subs=list(itertools.combinations(range(n),k))
            Pk=[[toZ(per([[E[r][u] for u in U] for r in R]),MOD) for U in subs] for R in subs]
            v=[Fm[sum(1<<u for u in U)]%MOD for U in subs]
            if not member(Pk,v): ok=False
    print(label,"degree-k parts in rowspace(Per_k):",ok,flush=True)
for (a,t) in []:
    E,c,T,B=cauchy_block(a,t,1,0); test(E,f"block a={a} t={t}")
# control: Cauchy-like without H
X=[0,2,4,6]; Y=[5,7,9,11]; c=Fr(1)
test([[c/(x-y) for y in Y] for x in X],"control no-H (c=1)")
X=[0,2,4,6]; Y=[5,7,9,11]; c=Fr(8)
test([[c/(x-y) for y in Y] for x in X],"H holds (c=8)")
