# Scalar-level (eigen-sign-pattern) check for m = 9 (Lucas, a = 2) and m = 3 with larger Q, over GF(27).
import itertools, math, sys
from gf import GF
from design import *
F=GF(3)
for (Q,m) in [(10,9),(12,9),(13,9),(12,3),(6,9)]:
    taus=list(range(1,Q+1))      # 13 <= 26 distinct nonzero elements of GF(27)
    best={}
    for shift in range(1,m+1):
        des=Design(F,Q,m,shift,taus)
        maxdeg,nbp,acc=des.scalar_checks()
        fr={r:float(des.bad_fraction(r)) for r in (2,3,4,8)}
        for r in fr: best[r]=min(best.get(r,1),fr[r])
        print(f"Q={Q} m={m} shift={shift} blocks={[ (b[0],b[-1]) for b in des.blocks]} deg={maxdeg} badpairs={nbp} accidentally_good={acc} frac(r=2,3,4,8)={[round(fr[r],4) for r in fr]}")
    print(f"== Q={Q} m={m}: best-shift error fraction by r: { {r:round(v,4) for r,v in best.items()} } ; 1/m={1/m:.4f}")
# solvability of (#): for every block of m consecutive integers (m=3,9,27 by linear algebra mod 3), and failure for m not a power of 3
def solvable(block,m):
    import numpy as np
    A=[[math.comb(w,i)%3 for i in range(m)]+[((-1)**w)%3] for w in block]
    A=np.array(A); rk=0
    rows,cols=A.shape
    for c in range(cols-1):
        piv=[i for i in range(rk,rows) if A[i,c]%3]
        if not piv: continue
        A[[rk,piv[0]]]=A[[piv[0],rk]]
        if A[rk,c]==2: A[rk]=(2*A[rk])%3
        for i in range(rows):
            if i!=rk and A[i,c]: A[i]=(A[i]-A[i,c]*A[rk])%3
        rk+=1
    return not any((not A[i,:-1].any()) and A[i,-1] for i in range(rows))
for m in (3,9,27,81):
    ok=all(solvable(list(range(lo,lo+m)),m) for lo in range(0,200))
    print(f"(#) solvable for all blocks [lo,lo+{m}), lo<200: {ok}")
for m in (2,4,5,6,7,8,10,12):
    bad=[lo for lo in range(0,60) if not solvable(list(range(lo,lo+m)),m)]
    print(f"m={m}: unsolvable full blocks start at {bad[:8]}{'...' if len(bad)>8 else ''}")
