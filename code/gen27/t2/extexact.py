# exact block-diagonal extension N0 -> N0+1 at the same q: E'_l = E_l (+) e_l, C'_l = C_l (+) c_l.
# Cross cells cost N0 iff M1 = sum c_l E_l, M2 = sum e_l C_l satisfy M1 M2 = I; cell (new,new) free iff <e,c> = 1.
# Search c in {-1,0,1}^q, M1 invertible, M2 = M1^{-1} must lie in span(C) with coefficient e, check <e,c> = 1 (else cost +1).
import sys, itertools, numpy as np, sympy
from lib import *
def extend(N0,E,C,maxc=1):
    q=len(E); best=None
    Cm=sympy.Matrix([[int(v) for v in c.ravel()] for c in C]).T   # N0^2 x q
    vals=range(-maxc,maxc+1)
    for cv in itertools.product(vals,repeat=q):
        if not any(cv): continue
        M1=sum(ci*e for ci,e in zip(cv,E))
        d=round(np.linalg.det(M1.astype(float)))
        if d==0: continue
        M1s=sympy.Matrix(M1.tolist()); M2=M1s.inv()
        if any(x.q!=1 for x in M2): continue   # want integer
        try:
            sol,params=Cm.gauss_jordan_solve(sympy.Matrix(list(M2)))
        except ValueError: continue
        if params.shape[0]: sol=sol.subs({t:0 for t in params})
        ev=[int(x) for x in sol]
        if any(not float(x).is_integer() for x in sol): continue
        ip=sum(a*b for a,b in zip(ev,cv))
        cost=0 if ip==1 else 1
        if best is None or cost<best[0]: best=(cost,ev,list(cv))
        if cost==0: break
    return best
def apply(N0,E,C,ev,cv):
    N=N0+1; E2=[];C2=[]
    for l in range(len(E)):
        a=np.zeros((N,N),dtype=int); a[:N0,:N0]=E[l]; a[N0,N0]=ev[l]; E2.append(a)
        b=np.zeros((N,N),dtype=int); b[:N0,:N0]=C[l]; b[N0,N0]=cv[l]; C2.append(b)
    return N,E2,C2
if __name__=='__main__':
    f=sys.argv[1]; N0,E,C=load(f); r=extend(N0,E,C)
    print(f, r)
    if r:
        N,E2,C2=apply(N0,E,C,r[1],r[2]); print("check p,q:",pq(N,E2,C2))
