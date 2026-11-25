# G25-P10: effect of a cheaper leaf contraction (sweep / cut-width) on App. D row G (P3 model).
import sys, math, time; sys.path.insert(0,'.')
from model2 import run, W_INC
from scipy.optimize import minimize
def optG(lam, D, leaffac):
    def mk(y):
        u0=math.exp(y[0]); rho=math.exp(y[1])
        us=[max(2,math.floor(u0*rho**d)) for d in range(D)]
        cs=[max(0,math.ceil(y[2]+d*y[3])) for d in range(D)]
        return us,cs
    f=lambda y: run(lam,*mk(y),w=W_INC,kappa='lin',pack='slack',Dmode='adapt',amode='slack',r=64,leaf='depth',leaffac=leaffac)
    best=None
    for g0 in (lam/27,lam/40):
        res=minimize(f,[math.log(g0),math.log(0.1),12.0,2.0],method='Nelder-Mead',options={'maxiter':600,'xatol':1e-4,'fatol':1e-9})
        if best is None or res.fun<best.fun: best=res
    return best.fun
def thr(tau, leaffac):
    lam=1e6
    for it in range(5):
        best=min(optG(lam,D,leaffac) for D in (1,2,3,4))
        new=(best*math.sqrt(lam)/tau)**2
        done=abs(new/lam-1)<0.003; lam=new
        if done: break
    return lam
t0=time.time()
for lf in (1.0,0.5,0.25,0.0):
    E=min(optG(8.97e5,D,lf) for D in (1,2,3,4))
    print(f"leaffac={lf}: E/n at lam=8.97e5: {E:.4f}; threshold vs 2^n: {thr(1.0,lf):.4g}  [{time.time()-t0:.0f}s]",flush=True)
