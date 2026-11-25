# G26-P8: recompute rows G/H thresholds with P3's full finite-size model (model2.run), bisection, densities 4.26 and 1.
import sys, math, time
sys.path.insert(0,'.')
import model2
from scipy.optimize import minimize
K=3
def set_alpha(alpha):
    model2.LOGN_OFF=math.log2(1+alpha)
    return K*alpha     # W_INC = km/n
def optG(lam, D, w, r=64, smith=False, kappa='lin', ret=False, starts=((1/27,0.1),(1/40,0.1),(1/60,0.2))):
    def mk(y):
        u0=math.exp(y[0]); rho=math.exp(y[1])
        us=[max(2,math.floor(u0*rho**d)) for d in range(D)]
        cs=[max(0,math.ceil(y[2]+d*y[3])) for d in range(D)]
        return us,cs
    f=lambda y: model2.run(lam,*mk(y),w=w,kappa=kappa,pack='slack',Dmode='adapt',amode='slack',r=r,leaf='depth',smith=smith)
    best=None
    for g0,gr in starts:
        res=minimize(f,[math.log(lam*g0),math.log(gr),12.0,2.0],method='Nelder-Mead',options={'maxiter':800,'xatol':1e-5,'fatol':1e-10})
        if best is None or res.fun<best.fun: best=res
    if ret: return best.fun, mk(best.x)
    return best.fun
def Ebest(lam,w,Ds=(2,3,4),**kw):
    return min((optG(lam,D,w,**kw),D) for D in Ds)
def thr(tau,w,lo,hi,**kw):
    # geometric bisection on E(lam)/n < tau
    while hi/lo>1.002:
        mid=math.sqrt(lo*hi)
        if Ebest(mid,w,**kw)[0]<tau: hi=mid
        else: lo=mid
    return hi
if __name__=='__main__':
    t0=time.time()
    which=sys.argv[1] if len(sys.argv)>1 else 'all'
    for alpha in [float(sys.argv[1])]:
        w=set_alpha(alpha)
        for nm,kw in [r for r in (("G",{}),("H",dict(smith=True))) if r[0] in sys.argv[2]]:
            for tau,tn in [t for t in ((1.0,'2^n'),(0.386,'PPSZ')) if t[1] in sys.argv[3].split(',')]:
                # bracket from closed form c/sqrt(lam)*w = tau
                c=74.0 if nm=='G' else 37.0
                g=(c*w/tau)**2
                lo,hi=g/1.4,g*1.4
                assert Ebest(lo,w,**kw)[0]>=tau and Ebest(hi,w,**kw)[0]<tau
                L=thr(tau,w,lo,hi,**kw)
                E,D=Ebest(L,w,**kw)
                print(f"alpha={alpha} row {nm} vs {tn}: lam*={L:.4g}  (E/n={E:.4f} D={D})  E*sqrt(lam)/(km)={E*math.sqrt(L)/w:.2f}  [{time.time()-t0:.0f}s]",flush=True)
