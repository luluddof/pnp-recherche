# Levenberg-Marquardt feasibility test for (N,p,q), real or complex.  [N] numerical only:
# "feasible" = residual < 1e-10 found; "not found" is NOT a proof of infeasibility.
import numpy as np, sys, time
from scipy.optimize import least_squares
def crossing(N):
    T = np.zeros((N,N,N,N))
    for x in range(N):
        for y in range(N): T[x,y,y,x] = 1.0
    return T
def model(z, N, p, q, cplx):
    n2 = N*N; m = (2*p+2*q)*n2
    if cplx: z = z[:m] + 1j*z[m:]
    F = z[:p*n2].reshape(p,N,N); Fp = z[p*n2:2*p*n2].reshape(p,N,N)
    G = z[2*p*n2:(2*p+q)*n2].reshape(q,N,N); Gp = z[(2*p+q)*n2:].reshape(q,N,N)
    return np.einsum('jxy,juv->xyuv', F, Fp) + np.einsum('kxu,kyv->xyuv', G, Gp)
def test(N, p, q, cplx, restarts, rng, tmax):
    T = crossing(N); n2 = N*N; m = (2*p+2*q)*n2; best = 1e9; t0 = time.time(); sol = None
    def fun(z):
        d = (model(z, N, p, q, cplx) - T).ravel()
        return np.concatenate([d.real, d.imag]) if cplx else d
    for r in range(restarts):
        z0 = rng.standard_normal(m*(2 if cplx else 1)) * rng.choice([0.3, 1.0])
        s = least_squares(fun, z0, method='trf' if fun(z0).size < z0.size else 'lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
        res = np.sqrt(2*s.cost)
        if res < best: best = res; sol = s.x
        if best < 1e-10 or time.time()-t0 > tmax: break
    return best, r+1, (np.abs(sol).max() if sol is not None else None)
if __name__ == "__main__":
    rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
    for spec in sys.argv[1].split(';'):
        N, p, q, f = spec.split(','); N, p, q = int(N), int(p), int(q)
        lbp = max(k*(N-k*q) for k in range(N+1)); lbq = max(k*(N-k*p) for k in range(N+1))
        best, r, mx = test(N, p, q, f == "c", 40, rng, 28)
        print(f"N={N} (p,q)=({p},{q}) field={f.upper()} bound_ok={p>=lbp and q>=lbq} best_resid={best:.2e} restarts={r} max|coef|={mx:.1f} ->", "FEASIBLE(num)" if best < 1e-10 else "not found", flush=True)
