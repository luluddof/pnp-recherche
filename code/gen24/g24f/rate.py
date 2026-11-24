# per-copy rate bounds for A1^n, xi = (alpha/2)^n, labels M = V' (or coset), via generating functions (rigorous Chernoff)
import math, numpy as np
from scipy.optimize import minimize
SQ8=2*math.sqrt(2)
K=np.arange(1,3000)
def lphi(b): return -np.log1p(-np.exp(-b*K)).sum()           # log 1/prod(1-e^{-bk})
S=np.arange(-12,14)
def wsec(s,shift=0.0): return (s+shift)**2-(s+shift)          # sector weight n^2-n (shift=1/2: coset alpha/2+L)
MSHIFT=0.0
def logZD(b):  # Dip profile (V' or M)
    w=wsec(S); return np.log(np.exp(-b*w).sum())+lphi(b)
def logZM(b):
    w=wsec(S,MSHIFT); w=w-w.min()*0; return np.log(np.exp(-b*w).sum())+lphi(b)
def leg(f,x,b0=1.0):
    r=minimize(lambda t: f(math.exp(t[0]))+math.exp(t[0])*x,[math.log(b0)],method='Nelder-Mead',options={'xatol':1e-7,'fatol':1e-10})
    return r.fun
def logBC(bx,by,lam=0.5):
    # sum_q R_q^{1-lam} R_{q-1}^lam,  R_q = sum_s e^{-by w(s)} e^{-bx wM(q+s)}  (times 1/phi phi)
    Q=np.arange(-24,26)
    ws=wsec(S)
    R=np.array([np.exp(-by*ws - bx*wsec(q+S,MSHIFT)).sum() for q in Q])
    Rm=np.array([np.exp(-by*ws - bx*wsec(q-1+S,MSHIFT)).sum() for q in Q])
    return math.log((R**(1-lam)*Rm**lam).sum())+lphi(bx)+lphi(by)
def flat(u,v,lam=0.5):
    f=lambda t: logBC(math.exp(t[0]),math.exp(t[1]),lam)+math.exp(t[0])*u+math.exp(t[1])*v
    best=None
    for t0 in [(0,0),(-1,-1),(1,1),(-1,1),(1,-1)]:
        r=minimize(f,list(t0),method='Nelder-Mead',options={'xatol':1e-7,'fatol':1e-10,'maxiter':2000})
        if best is None or r.fun<best: best=r.fun
    return best
if __name__=="__main__":
    import sys
    if len(sys.argv)>1: MSHIFT=float(sys.argv[1])
    sD=lambda x: leg(logZD,x); sM=lambda x: leg(logZM,x)
    grid=np.exp(np.linspace(math.log(0.05),math.log(4.0),40))
    best={'loose':(9,),'flat':(9,)}
    rows=[]
    for u in grid:
        for v in grid:
            cost=2*sM(u+v)
            loose=min(sM(u),sD(v),sM(2*u+2*v)/2)
            fl=flat(u,v)/2
            r1=cost/loose; r2=cost/min(loose,fl)
            if r1<best['loose'][0]: best['loose']=(r1,u,v)
            if r2<best['flat'][0]: best['flat']=(r2,u,v,loose,fl)
    print("MSHIFT",MSHIFT,"loose min 4rho=%.4f at u=%.3f v=%.3f"%best['loose'])
    print("with flattening: min 4rho=%.4f at u=%.3f v=%.3f (loose lam=%.4f flat lam=%.4f)"%best['flat'])
