# per-copy (tensor-power family) criterion, Legendre entropies as in G22e modscale.py
import sys, math, numpy as np
sys.path.insert(0,'/home/user/pnp-recherche/code/gen22/g22e')
from chars import mul, heis
from scaling import mk_phi
N=300
def entropy(d,den=1,mult0=1.0):
    phi,mean=mk_phi(d,den)
    bs=np.exp(np.linspace(math.log(40),math.log(0.02),2500))
    B=[b for b in bs if mean(b)[1]<1e-12]
    E=np.array([mean(b)[0] for b in B]); S=np.array([phi(b)+math.log(mult0)+b*e for b,e in zip(B,E)])
    o=np.argsort(E); E,S=E[o],S[o]
    s0=math.log(d[0]*mult0)
    return (lambda e: np.interp(e,np.concatenate([[0],E]),np.concatenate([[s0],S]))), E.max()
def lat(f,N=N):
    num=[0]*(N+1)
    for n in range(-40,41):
        e=f(n)
        if 0<=e<=N: num[e]+=1
    return mul(num,heis(N),N)
def crit(sCost, lam_bounds, emax=7.0, tag="", ngrid=150):
    """min over (u=e, v=h/2) of 2 sCost(u+v)/min(bounds(u,v)), bounds: list of functions (u,v)->rate"""
    us=np.concatenate([[0],np.exp(np.linspace(math.log(1e-3),math.log(emax/2),ngrid))]); bv=(9,0,0)
    for u in us:
        for v in us:
            if u+v>emax: continue
            lam=min(f(u,v) for f in lam_bounds)
            if lam>1e-9:
                r=2*sCost(u+v)/lam
                if r<bv[0]: bv=(r,u,v)
    print("%-62s min 4rho=%.4f at e=%.3g h/2=%.3g"%(tag,*bv),flush=True); return bv
if __name__=="__main__":
    Dip=lat(lambda n:n*n-n); Sym=lat(lambda n:n*n+abs(n)); Conf=lat(lambda n:n*n)
    sD,_=entropy(Dip); sS,_=entropy(Sym); sC,_=entropy(Conf); sDh,_=entropy(Dip,1,0.5)
    U=lambda u,v: sD(2*u+2*v)/2
    crit(sD,[lambda u,v:sD(u),lambda u,v:sD(v),U],tag="(d) loose: labels Dip, pairs Dip  [G22e 2.8157]")
    crit(sD,[lambda u,v:sD(u),lambda u,v:sS(v),U],tag="(a) channel 0 only: labels Dip, pairs Sym [G22e 3.030]")
    crit(sD,[lambda u,v:sS(u),lambda u,v:sD(v),U],tag="(a') channel alpha only: labels Sym, pairs Dip")
    crit(sD,[lambda u,v:sD(u),lambda u,v:sC(v),U],tag="(a'') channel 0, pairs by conformal lemma [G22e 2.916]")
    crit(sD,[lambda u,v:sD(u)-math.log(2),lambda u,v:sD(v),U],tag="(c) branch-2 product (labels Dip/2 per copy, pairs Dip)")
