"""min rank(P - u0 (x) u1 (x) u2) at N=2 (bond dimensions (1,1,1)); many restarts; real and complex."""
import sys, numpy as np, time
from hexlib import *
N = 2; d = 4; P = perm_P(N)
rs = [int(a) for a in sys.argv[1].split(',')]; starts = int(sys.argv[2])
t0 = time.time()
for r in rs:
    fg, unpack, n = make_obj(P, d, (1, 1, 1), r)
    rng = np.random.default_rng(100 + r); vals = []; best = (9e9, None)
    for s in range(starts):
        sc = [0.3, 0.6, 1.0][s % 3]
        x0 = rng.normal(size=n) * sc
        if s % 4 == 1:   # near +-identity / local unitaries
            x0 = np.concatenate([(np.eye(4) * rng.choice([-1, 1]) + 0.3 * rng.normal(size=(4, 4))).ravel() for _ in range(3)])
        res = minimize(fg, x0, jac=True, method='L-BFGS-B', options=dict(maxiter=3000, maxfun=6000, ftol=1e-16, gtol=1e-11))
        vals.append(res.fun)
        if res.fun < best[0]: best = (res.fun, res.x.copy())
        if res.fun < 1e-10: break
    vals = np.sort(vals)
    print(f"r={r}: best residual^2={best[0]:.3e}; lowest minima {np.round(vals[:6], 4)}  starts={len(vals)} [{time.time()-t0:.0f}s]", flush=True)
    if best[0] < 1e-8:
        B0, B1, B2 = unpack(best[1]); np.save(f"prodB_r{r}.npy", best[1])
        s = np.linalg.svd(P - build(B0, B1, B2), compute_uv=False)
        print("   exact-rank check: sigma_r, sigma_{r+1} =", s[r-1], s[r], " |x|=", np.linalg.norm(best[1]))
