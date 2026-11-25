"""Precise alpha(nu) for symmetric frontiers, nu in [1/2,1): excess(nu) = h(nu) - nu = max_t (1-nu) f(1-t) - nu t,
condition for alpha: excess < nu * expm1(log1p(r^q)/q), r = (1-nu)/nu, q = 1/(1-alpha)."""
import numpy as np, sys, math
sys.path.insert(0, '/tmp/claude-0/-home-user-pnp-recherche/2d12e629-e6d2-542b-bfc0-8b766f57383f/scratchpad/gen29/u1')
from sim import families
t = np.unique(np.concatenate([np.geomspace(1e-300, 1, 400001), np.linspace(0, 1, 200001)]))
def alpha_nu(exc, nu):
    r = (1-nu)/nu
    def ok(a):
        if a >= 1: return exc <= 0
        q = 1/(1-a); return exc < nu*math.expm1(math.log1p(r**q)/q)
    if not ok(0.0): return 0.0
    lo, hi = 0.0, 1.0
    for _ in range(80):
        m = 0.5*(lo+hi)
        if ok(m): lo = m
        else: hi = m
    return lo
fam = families()
fam['circle'] = lambda u: np.sqrt(np.clip((1-u)*(1+u), 0, None))
nus = np.unique(np.concatenate([np.linspace(0.5, 0.999, 500), 1 - np.geomspace(1e-12, 1e-3, 200)]))
for nm, f in fam.items():
    best = (0, 0); rows = []
    for nu in nus:
        u = 1 - t
        fv = f(u) if nm != 'circle' else np.sqrt(np.clip(t*(2-t), 0, None))
        exc = max(0.0, float(np.max((1-nu)*fv - nu*t)))
        a = alpha_nu(exc, nu); rows.append((nu, a))
        if a > best[0]: best = (a, nu)
    sel = [r for r in rows if abs(r[0]-0.5) < 1e-9 or abs(r[0]-0.8) < 1e-3 or abs(r[0]-0.95)<1e-3 or r[0] > 1-2e-12 or abs(r[0]-(1-1e-6))<1e-7]
    print('%-14s alpha*=%.4f (nu=%.6g)  ' % (nm, best[0], best[1]) + ' '.join('a(%.6g)=%.4f' % r for r in sel))
