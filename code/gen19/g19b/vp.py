# Variable projection: minimise over C = sum_{l<=q} G_l (x) G'_l the tail  t(C) = sum_{i>p} sigma_i(F - C)^2  (row flattening).
# t = 0 <=> (p,q) feasible.  L-BFGS, random restarts.  [N] numerical only.
import numpy as np, sys, time
from scipy.optimize import minimize
def swap(N):
    F = np.zeros((N*N, N*N))
    for x in range(N):
        for y in range(N): F[x*N+y, y*N+x] = 1.0
    return F
def fg(z, N, p, q, F):
    G = z[:q*N*N].reshape(q,N,N); Gp = z[q*N*N:].reshape(q,N,N)
    C = np.einsum('kxu,kyv->xyuv', G, Gp).reshape(N*N, N*N)
    U, s, Vt = np.linalg.svd(F - C)
    Et = (U[:, p:] * s[p:]) @ Vt[p:]
    Et4 = Et.reshape(N,N,N,N)
    gG = -2*np.einsum('xyuv,kyv->kxu', Et4, Gp); gGp = -2*np.einsum('xyuv,kxu->kyv', Et4, G)
    return float((s[p:]**2).sum()), np.concatenate([gG.ravel(), gGp.ravel()])
def run(N, p, q, tmax, seed, maxit=4000):
    F = swap(N); rng = np.random.default_rng(seed); t0 = time.time(); hist = []; best = (1e9, None)
    while time.time() - t0 < tmax:
        z0 = rng.standard_normal(2*q*N*N) * rng.choice([0.3, 0.6, 1.0])
        r = minimize(fg, z0, args=(N,p,q,F), jac=True, method='L-BFGS-B', options=dict(maxiter=maxit, maxfun=maxit*2, ftol=1e-30, gtol=1e-14))
        hist.append((r.fun, float(np.abs(r.x).max())))
        if r.fun < best[0]: best = (r.fun, r.x.copy())
        if best[0] < 1e-22: break
    return best, hist
if __name__ == "__main__":
    tmax = float(sys.argv[2]); maxit = int(sys.argv[3]) if len(sys.argv) > 3 else 4000
    for spec in sys.argv[1].split(';'):
        N, p, q = map(int, spec.split(','))
        (t, z), hist = run(N, p, q, tmax, 1000*N + 10*p + q, maxit)
        hs = sorted(hist)[:6]
        tag = "FEASIBLE(num)" if t < 1e-20 else ("border?" if t < 1e-3 else "not found")
        print(f"N={N} (p,q)=({p},{q}) best tail={t:.3e} max|G|={np.abs(z).max():.1f} restarts={len(hist)} lowest={[(float('%.2e'%a), float('%.0f'%b)) for a,b in hs]} -> {tag}", flush=True)
        np.save(f"vp_{N}_{p}_{q}.npy", z)
