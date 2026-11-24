# vp.py with warm starts: initial C = (solution file for N0 < N, zero-padded to N, extra terms random) + noise.  [N]
import numpy as np, sys, time
from scipy.optimize import minimize
from vp import swap, fg
N, p, q = map(int, sys.argv[1].split(',')); init = sys.argv[2]; N0, q0 = map(int, sys.argv[3].split(','))
tmax = float(sys.argv[4]); maxit = int(sys.argv[5]); seed = int(sys.argv[6]) if len(sys.argv) > 6 else 0
z = np.load(init); G0 = z[:q0*N0*N0].reshape(q0, N0, N0); H0 = z[q0*N0*N0:].reshape(q0, N0, N0)
F = swap(N); rng = np.random.default_rng(seed); t0 = time.time(); best = (1e9, None); hist = []
while time.time() - t0 < tmax:
    G = rng.standard_normal((q, N, N)) * 0.3; H = rng.standard_normal((q, N, N)) * 0.3
    k = min(q, q0); G[:k, :N0, :N0] = G0[:k] + 0.05*rng.standard_normal((k, N0, N0)); H[:k, :N0, :N0] = H0[:k] + 0.05*rng.standard_normal((k, N0, N0))
    r = minimize(fg, np.concatenate([G.ravel(), H.ravel()]), args=(N, p, q, F), jac=True, method='L-BFGS-B', options=dict(maxiter=maxit, maxfun=2*maxit, ftol=1e-30, gtol=1e-14))
    hist.append((float('%.3e' % r.fun), float('%.1f' % np.abs(r.x).max()))); print("restart", len(hist), hist[-1], f"{time.time()-t0:.0f}s", flush=True)
    if r.fun < best[0]: best = (r.fun, r.x.copy())
    if best[0] < 1e-22: break
t = best[0]; tag = "FEASIBLE(num)" if t < 1e-20 else ("border?" if t < 1e-3 else "not found")
print(f"N={N} (p,q)=({p},{q}) warm start from {init}: best tail={t:.3e} restarts={len(hist)} lowest={sorted(hist)[:6]} -> {tag}", flush=True)
np.save(f"vpi_{N}_{p}_{q}.npy", best[1])
