# Exact algebraic check of P8 Theorem U identity: N = Gamma * (+)_k Per_k(E), for E_xy = c/(X_x - Y_y), arbitrary rationals.
import itertools, random, sys
from fractions import Fraction as Fr
from lib12 import *
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 1)
def mobius_cols(M,n):
    N=1<<n; R=[row[:] for row in M]
    for j in range(n):
        for row in R:
            for X in range(N):
                if X>>j&1: row[X]-=row[X^(1<<j)]
    return R
def gamma(x,B,c,X):
    p=Fr(1)
    for xp in B:
        if xp!=x: p*=c/(X[xp]-X[x])
    return p
def Gam(S,R,c,X):
    Sl=[x for x in S if x not in R]; Rl=list(R); tot=Fr(0)
    for phi in itertools.product(Rl,repeat=len(Sl)):
        pre={r:[r] for r in Rl}
        for s,r in zip(Sl,phi): pre[r].append(s)
        p=Fr(1)
        for r in Rl: p*=gamma(r,pre[r],c,X)
        tot+=p
    return tot
bad=0
for trial in range(12):
    n=random.choice([2,3,4])
    X=[Fr(random.randint(-50,50),random.randint(1,9)) for _ in range(n)]
    while len(set(X))<n: X=[Fr(random.randint(-50,50),random.randint(1,9)) for _ in range(n)]
    Y=[Fr(random.randint(-50,50),random.randint(1,9))+Fr(1,997) for _ in range(n)]
    c=Fr(random.randint(1,30),random.randint(1,7))
    E=[[c/(X[x]-Y[y]) for y in range(n)] for x in range(n)]
    M=Mmat(E); Nm=mobius_cols(M,n)
    for S in range(1<<n):
        Ss=[x for x in range(n) if S>>x&1]
        for U in range(1<<n):
            Us=[y for y in range(n) if U>>y&1]
            rhs=Fr(0)
            for R in itertools.combinations(Ss,len(Us)):
                rhs+=Gam(Ss,R,c,X)*per([[E[r][y] for y in Us] for r in R])
            if rhs!=Nm[S][U]: bad+=1
    print("trial",trial,"n",n,"bad so far",bad,flush=True)
print("TOTAL mismatches",bad)
