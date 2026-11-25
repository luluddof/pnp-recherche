# Smith(M(E)) vs union_k Smith(Per_k E), and (A), and C1'-right, on real crossing blocks.
import sys, itertools, time
from fractions import Fraction as Fr
from lib12 import *
NB=int(sys.argv[2]) if len(sys.argv)>2 else 160
MOD=1<<NB
def run(a,t,zeta,b):
    E,c,T,B=cauchy_block(a,t,zeta,b); n=len(E)
    # integrality of gamma pairs
    gmin=min(v2q(c/(xp-x)) for x in B for xp in B if x!=xp) if n>1 else None
    M=Mmat(E)
    sm=smith_vals([[toZ(z,MOD) for z in r] for r in M],NB)
    un=[]; okA=True; sqf=[]
    sE=smith_vals([[toZ(z,MOD) for z in r] for r in E],NB)
    for k in range(n+1):
        if k==0: sp=[0]; sc=[0]
        else:
            P=compound(E,k,per); sp=smith_vals([[toZ(z,MOD) for z in r] for r in P],NB)
            sc=sorted(sum(cmb) for cmb in itertools.combinations(sE,k))
        un+=sp; sqf+=sc
        if any(p>q for p,q in zip(sp,sc)): okA=False
    un.sort(); sqf.sort()
    U = (sm==un)
    C1r = all(p<=q for p,q in zip(sm,sqf))
    print(f"a={a} t={t} zeta={zeta} b={b} n={n} min v2 gamma-factor={gmin} U={U} A={okA} C1'right={C1r} maxSmith={max(sm)}",flush=True)
    if not U: print("   M:",sm[:20],"\n   un:",un[:20])
    return U,okA,C1r
a=int(sys.argv[1]); L=1<<a
res=[]
zs=[1,-3,Fr(1,3),2,-12,Fr(5,7),8] if a<4 else [1,Fr(-5,3),2]
for t in sorted(set([t for t in range(-L,L) if t!=0 and v2(t)<a])):
    import os
    TL=os.environ.get("TL")
    if a==4 and t not in (eval(TL) if TL else (1,-3,2,-6,4,12,-8)): continue
    for zeta in zs:
        for b in (0,1):
            if a==4 and b==1 and zeta!=1: continue
            res.append(((a,t,zeta,b),run(a,t,zeta,b)))
