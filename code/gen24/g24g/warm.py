# [N] warm start: embed an exact N0 solution (json) into N = N0 + d (zero padding + random extension), run real vp for (p,q).
import numpy as np, sys, json, time
from scipy.optimize import minimize
sys.argv += [] 
src = json.load(open(sys.argv[1])); N = int(sys.argv[2]); p = int(sys.argv[3]); R = int(sys.argv[4]); amp = float(sys.argv[5]); MI = int(sys.argv[6])
E0 = np.array(src["E"], float); C0 = np.array(src["C"], float); q, N0 = E0.shape[0], E0.shape[1]
exec(open('vpreal.py').read().split('N,p,q,R,seed=')[0])   # swap, fg
F = swap(N); rng = np.random.default_rng(int(time.time()) % 1000); best = 1e9
for r in range(R):
    G = amp*rng.standard_normal((q, N, N)); H = amp*rng.standard_normal((q, N, N))
    G[:, :N0, :N0] = E0; H[:, :N0, :N0] = C0
    z0 = np.concatenate([G.ravel(), H.ravel()])
    res = minimize(lambda z: fg(z, F, N, p, q), z0, jac=True, method='L-BFGS-B', options=dict(maxiter=MI, ftol=1e-16, gtol=1e-12, maxcor=30))
    print(f"warm N={N} p={p} q={q} r{r} d2={res.fun:.3e} coef={np.abs(res.x).max():.2f} nit={res.nit}", flush=True)
    if res.fun < best: best = res.fun; np.save(f"warm_{N}_{p}_{q}.npy", res.x)
print("BEST", best)
