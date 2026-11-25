import sys, math; sys.path.insert(0,'.')
from model2 import *
from scipy.optimize import minimize
def optG(lam, D, r=64, smith=False, kappa='lin', verbose=False):
    def f(y):
        u0=math.exp(y[0]); rho=math.exp(y[1]); c=y[2]
        us=[max(2,math.floor(u0*rho**d)) for d in range(D)]
        cs=[max(0,math.ceil(c+ d*y[3])) for d in range(D)]
        return run(lam, us, cs, w=W_INC, kappa=kappa, pack='slack', Dmode='adapt', amode='slack', r=r, leaf='depth', smith=smith)
    best=None
    for g0 in (lam/25, lam/50, lam/100):
        for gr in (0.07,0.15,0.3):
            res=minimize(f,[math.log(g0),math.log(gr),12.0,2.0],method='Nelder-Mead',options={'maxiter':4000,'xatol':1e-5,'fatol':1e-10})
            if best is None or res.fun<best.fun: best=res
    if verbose:
        y=best.x; u0=math.exp(y[0]); rho=math.exp(y[1])
        us=[max(2,math.floor(u0*rho**d)) for d in range(D)]; cs=[max(0,math.ceil(y[2]+d*y[3])) for d in range(D)]
        print(f"   lam={lam:.3g} D={D} r={r}: E/n={best.fun:.5f}  E*sqrt(lam)/M0={best.fun*math.sqrt(lam)/W_INC:.2f}  us={us} cs={cs}")
        run(lam, us, cs, w=W_INC, kappa=kappa, pack='slack', Dmode='adapt', amode='slack', r=r, leaf='depth', smith=smith, verbose=True)
    return best.fun
if __name__=='__main__':
    for lam in (1e6,):
        for r in (32,64,128):
            for D in (2,3,4):
                optG(lam,D,r=r,verbose=True)
