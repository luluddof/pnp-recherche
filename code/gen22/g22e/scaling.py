# tensor-power scaling families V^{(x)n}, n -> infinity, h = eps*n : per-copy profile s(eps) = inf_b [phi(b)+b*eps]
# criterion (*) per copy: r(eps) = 2 s(2eps)/min(s(eps), s(4eps)/2) >= 2 sqrt2 ;  psi(b) = d/db[b phi(b)] sign scan
import math, numpy as np
from chars import *
SQ2=math.sqrt(2)
def mk_phi(d,den):
    w=np.array([i/den for i in range(len(d))]); ld=np.array([math.log(x) if x>0 else -np.inf for x in d])
    def phi(b):
        t=ld-b*w; m=t.max(); return m+math.log(np.exp(t-m).sum())
    def mean(b):
        t=ld-b*w; m=t.max(); p=np.exp(t-m); p/=p.sum(); return (p*w).sum(), p[-int(len(d)*0.2):].sum()
    return phi,mean
def run(name,d,den=1,mult0=1.0):
    phi0,mean=mk_phi(d,den)
    phi=lambda b: phi0(b)+math.log(mult0)
    # valid b: tail mass in last 20% of truncation < 1e-12
    bs=np.exp(np.linspace(math.log(30),math.log(0.02),3000))
    good=[b for b in bs if mean(b)[1]<1e-12]
    bmin=min(good)
    B=np.array(sorted(good))
    P=np.array([phi(b) for b in B])
    # psi = phi + b phi'
    E=np.array([mean(b)[0] for b in B])   # eps(b) = -phi'(b)
    psi=P-B*E
    sgn=np.sign(psi); changes=[(B[i],sgn[i],sgn[i+1]) for i in range(len(B)-1) if sgn[i]!=sgn[i+1]]
    # s(eps) on the valid range: s = phi(b)+b eps at eps=E(b)
    S=P+B*E
    order=np.argsort(E); E1=E[order]; S1=S[order]
    sf=lambda e: np.interp(e,E1,S1)
    emax=E1.max(); eps=np.exp(np.linspace(math.log(E1.min()),math.log(emax/4),400))
    r=[2*sf(2*e)/min(sf(e),sf(4*e)/2) for e in eps]
    i=int(np.argmin(r))
    print("%-34s b_min %.3f eps_max %.1f | psi sign changes %s | min r = %.4f at eps=%.3g (r at eps_max/4 %.4f)"%(
        name,bmin,emax,[(round(float(b),3),int(a),int(c)) for b,a,c in changes][:4],r[i],eps[i],r[-1]))
    return eps,r
if __name__=="__main__":
    N=400
    run("Heisenberg",heis(N))
    run("NS fermion",NSferm(2*N),2)
    run("Ramond fermion (ground sqrt2)",ramond(N),1,math.sqrt(2))
    run("E8 lattice",E8(N))
    run("Leech",leech(300))
    run("Moonshine",J_dims(300))
    run("Yang-Lee vac",RR(N,{2,3}))
    run("Yang-Lee module",RR(N,{1,4}))
    run("M(2,7) vac",RR(N,{3,4},7))
    run("M(2,9) vac",RR(N,{4,5},9))
    run("symplectic fermions",symplf(N))
    run("Heisenberg Z2-twisted",prod_series([(k,1) for k in range(1,2*N+1,2)],2*N),2)
