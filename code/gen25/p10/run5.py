import sys, math, time; sys.path.insert(0,'.')
from model2 import run, W_INC
from scipy.optimize import minimize
def optG(lam, D, r=64, smith=False, kappa='lin', ret=False):
    def mk(y):
        u0=math.exp(y[0]); rho=math.exp(y[1])
        us=[max(2,math.floor(u0*rho**d)) for d in range(D)]
        cs=[max(0,math.ceil(y[2]+d*y[3])) for d in range(D)]
        return us,cs
    f=lambda y: run(lam,*mk(y),w=W_INC,kappa=kappa,pack='slack',Dmode='adapt',amode='slack',r=r,leaf='depth',smith=smith)
    best=None
    for g0 in (lam/27,lam/40):
        res=minimize(f,[math.log(g0),math.log(0.1),12.0,2.0],method='Nelder-Mead',options={'maxiter':600,'xatol':1e-4,'fatol':1e-9})
        if best is None or res.fun<best.fun: best=res
    if ret: return best.fun, mk(best.x)
    return best.fun
def thr(tau, **kw):
    lam=1e6
    for it in range(4):
        best=min(optG(lam,D,**kw) for D in (2,3))
        new=(best*math.sqrt(lam)/tau)**2
        done=abs(new/lam-1)<0.005; lam=new
        if done: break
    return lam
if __name__=='__main__':
    t0=time.time()
    rows=[("G  + a-slack, u_d∝M_d^{2/3}, depth 2-3, r=64 [rigorous]",{}),
          ("H  + Smith modes (nu/4) [conditional]",dict(smith=True)),
          ("I  + kappa=2 floor (black-box recovery) [open]",dict(smith=True,kappa=2.0))]
    for nm,kw in rows:
        a=thr(1.0,**kw); b=thr(0.386,**kw)
        c=optG(a,3,**kw)*math.sqrt(a)/W_INC
        print(f"{nm:60s} vs 2^n: {a:.3g}  vs PPSZ: {b:.3g}  E*sqrt(lam)/(km) @lam*: {c:.2f}   [{time.time()-t0:.0f}s]", flush=True)
    print("constant E*sqrt(lam)/(km), row G, r=64 best of D=2,3:")
    for lam in (1e5,1e6,1e7,1e8,1e9,1e10):
        v=min(optG(lam,D) for D in (2,3,4))
        print(f"  lam={lam:.0e}: {v*math.sqrt(lam)/W_INC:.2f}", flush=True)
    for r in (32,128):
        v=min(optG(1e6,D,r=r) for D in (2,3,4)); print(f"  r={r} lam=1e6: {v*1e3/W_INC:.2f}", flush=True)
