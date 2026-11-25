"""alpha* = inf{a : h(nu) >= ||(nu,1-nu)||_{1/(1-a)} for all nu}; also alpha(nu) per nu. Independent code."""
import numpy as np, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-pnp-recherche/2d12e629-e6d2-542b-bfc0-8b766f57383f/scratchpad/gen29/u1')
from sim import families

u = np.concatenate([np.linspace(0, 1, 200001), 1 - np.geomspace(1e-12, 1e-3, 20000)])
u = np.unique(np.clip(u, 0, 1))
nus = np.unique(np.concatenate([np.linspace(0.001, 0.999, 999), 1 - np.geomspace(1e-9, 1e-3, 300)]))

def qnorm(nu, q):
    m = np.maximum(nu, 1-nu); r = np.minimum(nu, 1-nu)/m
    return m*np.exp(np.log1p(r**q)/q)

def alpha_nu(hn, nu):
    # largest alpha with h(nu) < ||.||_{1/(1-alpha)}  (||.||_q decreasing in q)
    if hn >= qnorm(nu, 1.0) - 1e-15: return 0.0
    lo, hi = 0.0, 1.0
    for _ in range(60):
        a = 0.5*(lo+hi); q = 1/(1-a) if a < 1 else np.inf
        if a < 1 and hn < qnorm(nu, q): lo = a
        else: hi = a
    return lo

fam = families()
for nm, f in fam.items():
    v = f(u)
    an = []
    for nu in nus:
        hn = np.max(nu*u + (1-nu)*v)
        an.append(alpha_nu(hn, nu))
    an = np.array(an); k = int(np.argmax(an))
    i_half = int(np.argmin(abs(nus-0.5)))
    print('%-14s alpha*=%.4f at nu=%.6f | alpha(1/2)=%.4f (log2 1/rho) | alpha(nu=1-1e-9)=%.4f' % (
        nm, an[k], nus[k], an[i_half], an[-1]))
