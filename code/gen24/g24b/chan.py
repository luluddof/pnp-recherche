# single-channel bounds in the tensor-power family A1^n, xi=(alpha/2)^n, channel gamma with a fraction pi of alpha's
#   rigorous (avg) version: R_gamma <= 2 D_V^{(gamma/2)}(H/2) ~ Conf^{1-pi} Dip^{pi};  Lab_gamma <= D_M^{(xi-gamma/2)}(E0) ~ Dip^{1-pi} Conf^{pi}
#   exact-per-copy (max) version (heuristic for n copies): Sym in place of Conf
import math, numpy as np
from fam import *
from scaling import mk_phi
Dip=lat(lambda n:n*n-n); Sym=lat(lambda n:n*n+abs(n)); Conf=lat(lambda n:n*n)
bs=np.exp(np.linspace(math.log(30),math.log(0.03),1500))
def phis(d):
    phi,mean=mk_phi(d,1); return np.array([phi(b) for b in bs])
PD,PS,PC=phis(Dip),phis(Sym),phis(Conf)
def rate(P,x):  # inf_b [P(b)+b x] (cumulative count rate)
    return np.min(P+bs*x)
sD,_=entropy(Dip)
def run(Pother,tag):
    best=(9,)
    for pi in np.linspace(0,1,21):
        Pr=(1-pi)*Pother+pi*PD; Pl=(1-pi)*PD+pi*Pother
        us=np.exp(np.linspace(math.log(0.02),math.log(3.4),70)); bv=(9,0,0)
        for u in us:
            for v in us:
                lam=min(rate(Pr,v),rate(Pl,u),sD(u),sD(2*u+2*v)/2)
                r=2*sD(u+v)/lam
                if r<bv[0]: bv=(r,u,v)
        if bv[0]<best[0]: best=(bv[0],pi,bv[1],bv[2])
        print("  %s pi=%.2f: min 4rho=%.4f at e=%.3g h/2=%.3g"%(tag,pi,*bv),flush=True)
    print(tag,"BEST",best)
run(PC,"rigorous(Conf)")
run(PS,"max-version(Sym)")
