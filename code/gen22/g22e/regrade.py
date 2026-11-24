# A1 lattice VOA regraded by L0' = L0 - (1/2) alpha(0)  (sector n: weight n^2 - n + level).  Profile of V' = that of sl2_1 L(1/2) (has a dip).
# True pairing count at new weight h (pairs sector n with -n, a(2h-1)b in V'_0): <= sum_n dim(old weight h-|n|... ) = profile of Z_sym = sum q^{n^2+|n|}/eta
import math, numpy as np
from chars import *
from modscale import entropy
N=300
def lat(f):
    num=[0]*(N+1)
    for n in range(-40,41):
        e=f(n)
        if 0<=e<=N: num[e]+=1
    return mul(num,heis(N),N)
Zp=lat(lambda n:n*n-n); Zs=lat(lambda n:n*n+abs(n)); Z0=lat(lambda n:n*n)
sP,_=entropy(Zp); sS,_=entropy(Zs); s0,_=entropy(Z0)
def best(sM,sV,emax=7.0,tag=""):
    us=np.concatenate([[0],np.exp(np.linspace(math.log(1e-3),math.log(emax/2),150))]); bv=(9,0,0)
    for u in us:
        for v in us:
            if u+v>emax: continue
            lam=min(sM(u),sV(v),sM(2*u+2*v)/2)
            if lam>1e-9:
                r=2*sM(u+v)/lam
                if r<bv[0]: bv=(r,u,v)
    print("%-55s min 4rho=%.4f at u=%.3g v=%.3g"%(tag,*bv))
best(sP,s0,tag="regraded cost, anchor in CONFORMAL grading")
best(sP,sP,tag="regraded, LOOSE anchor (N<=dim V'_h)")
best(sP,sS,tag="regraded, TRUE pairing anchor")
best(s0,s0,tag="A1 lattice, conformal grading")
