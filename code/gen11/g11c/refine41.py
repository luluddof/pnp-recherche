# Exact-arithmetic refinement of the (4,1) profile-(2,3,3,3) solution: residuals evaluated EXACTLY (Fractions),
# Gauss-Newton steps in floating point.  tau = 1, 17/10, 29/10, 43/10.
import numpy as np
from fractions import Fraction as Fr
from feas import Prob
l,d=4,1; tauQ=[Fr(1),Fr(17,10),Fr(29,10),Fr(43,10)]; prof=[2,3,3,3]
P=Prob(l,d,[float(t) for t in tauQ],prof); n=1<<l
p=[Fr(float(v)) for v in np.load("sol41_2333.npy")]
cons=[]   # (k1,x,y,k2): (E_k1-X)/tau_k1 - (E_k2-X)/tau_k2 = 0
for x in range(n):
    for y in range(n):
        D=[k for k in range(l) if ((x^y)>>k)&1]
        for i in range(len(D)-1): cons.append((D[i],D[i+1],x,y))
def evalsQ(p):
    e=[None]*len(P.elist)
    for idx,(k,side,i,j) in enumerate(P.elist):
        r=prof[k]; oU=P.off[(k,side,'U')]; oV=P.off[(k,side,'V')]
        e[idx]=sum(p[oU+i*r+q]*p[oV+j*r+q] for q in range(r))
    return e
def resQ(p):
    e=evalsQ(p); out=[]
    for k1,k2,x,y in cons:
        X=1 if y==(n-1)^x else 0
        out.append((e[P.eidx[(k1,x,y)]]-X)/tauQ[k1]-(e[P.eidx[(k2,x,y)]]-X)/tauQ[k2])
    return out
C=np.zeros((len(cons),len(P.elist)))
for i,(k1,k2,x,y) in enumerate(cons):
    C[i,P.eidx[(k1,x,y)]]=1/float(tauQ[k1]); C[i,P.eidx[(k2,x,y)]]=-1/float(tauQ[k2])
def jacE(pf):
    W=P.W; P.W=np.eye(len(P.elist)); J=P.jac(pf); P.W=W; return J
for it in range(6):
    r=resQ(p); nr=max(abs(v) for v in r)
    print("iter",it,"exact max|residual| = %.3e"%float(nr),flush=True)
    if nr==0: break
    pf=np.array([float(v) for v in p]); J=C@jacE(pf)
    dp=np.linalg.lstsq(J,-np.array([float(v) for v in r]),rcond=1e-9)[0]
    p=[a+Fr(float(b)) for a,b in zip(p,dp)]
