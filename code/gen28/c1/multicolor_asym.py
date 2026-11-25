# h -> infinity limit of the counted multi-colour App-B curve, m = kappa*h.
# log p_m(x h) ~ h G(x),  G(x) = min_t [kappa L(t) + t x],  L(t) = sum_n -log(1-e^{-n t}).
import numpy as np
from scipy.optimize import minimize_scalar
def L(t):
    n=np.arange(1,4000); return float(np.sum(-np.log1p(-np.exp(-n*t))))
def G(x,k):
    if x<=0: return 0.0
    r=minimize_scalar(lambda lt: k*L(np.exp(lt))+np.exp(lt)*x, bounds=(-12,6), method='bounded', options={'xatol':1e-10})
    return r.fun
for k in [1e-4,1e-3,0.01,0.05,0.2,1,5]:
    g1=G(1,k); end=G(4,k)/g1
    # curve sigma in (-2,2): theta=G(2+s)/g1, T=G(2-s)/g1 ; jump = min T with theta<=2 ; balanced
    ss=np.linspace(-1.999,1.999,801)
    th=np.array([G(2+s,k)/g1 for s in ss]); T=np.array([G(2-s,k)/g1 for s in ss])
    Tstar=T[th<=2].min(); bal=np.min(np.maximum(th,T))
    # local exponent of cheap side near 0: G(x)/g1 at x=1e-3,1e-2
    loc=[G(x,k)/g1 for x in (1e-3,1e-2,1e-1)]
    print(f"kappa={k:g}: theta_end={end:.4f} (2 = corner reachable)  jump T*={Tstar:.4f}  balanced={bal:.4f}  C-side eta(x=1e-3,1e-2,1e-1)={loc[0]:.4f},{loc[1]:.4f},{loc[2]:.4f}")
