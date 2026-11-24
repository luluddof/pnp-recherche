# driver around G24c's vp.fg: per-restart printing, optional warm start from a smaller-N vp point (zero padding + new terms)
import numpy as np, sys, time
from scipy.optimize import minimize
from vp import fg, swap
N, p, q, R, seed, maxit = map(int, sys.argv[1:7]); warm = sys.argv[7] if len(sys.argv) > 7 else None
if q > p: p, q = q, p
F = swap(N).astype(complex); rng = np.random.default_rng(seed); s = N*N; best = 1e9
for r in range(R):
    if warm:
        N0, q0 = map(int, sys.argv[8:10]); z0 = np.load(warm); n0 = len(z0)//2; w0 = z0[:n0] + 1j*z0[n0:]
        G0 = w0[:q0*N0*N0].reshape(q0, N0, N0); H0 = w0[q0*N0*N0:].reshape(q0, N0, N0)
        G = np.zeros((q, N, N), complex); H = np.zeros((q, N, N), complex)
        perm = rng.permutation(N) if r % 2 else np.arange(N)     # also try random placements of the old block
        m = min(q, q0); G[:m][:, :N0, :N0] = G0[:m]; H[:m][:, :N0, :N0] = H0[:m]
        G = G[:, perm][:, :, perm]; H = H[:, perm][:, :, perm]
        sc = 0.3
        G += sc*(rng.standard_normal(G.shape) + 1j*rng.standard_normal(G.shape))*np.where(np.abs(G) > 0, 0.1, 1)
        H += sc*(rng.standard_normal(H.shape) + 1j*rng.standard_normal(H.shape))*np.where(np.abs(H) > 0, 0.1, 1)
        w = np.concatenate([G.ravel(), H.ravel()]); z0 = np.concatenate([w.real, w.imag])
    else:
        z0 = rng.standard_normal(4*q*s) * rng.choice([0.3, 0.7, 1.0, 1.5])
    t0 = time.time()
    res = minimize(lambda z: fg(z, F, N, p, q)[:2], z0, jac=True, method='L-BFGS-B',
                   options=dict(maxiter=maxit, maxfun=2*maxit, ftol=1e-16, gtol=1e-12, maxcor=30))
    f, g, S = fg(res.x, F, N, p, q)
    if f < best: best = f; np.save(f"b2_{N}_{p}_{q}_s{seed}{'w' if warm else ''}.npy", res.x)
    print(f"N={N} p={p} q={q} {'warm' if warm else 'cold'} r{r} d2={f:.6g} coef={np.abs(res.x).max():.2f} nit={res.nit} t={time.time()-t0:.0f}s sv={[round(float(t),4) for t in S[max(0,p-2):p+3]]}", flush=True)
