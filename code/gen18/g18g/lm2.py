# Same test as lm.py with an analytic Jacobian (fast).  [N] numerical only.
import numpy as np, sys, time
from scipy.optimize import least_squares
def crossing(N):
    T = np.zeros((N,N,N,N))
    for x in range(N):
        for y in range(N): T[x,y,y,x] = 1.0
    return T
def unpack(z, N, p, q, cplx):
    n2 = N*N; m = (2*p+2*q)*n2
    if cplx: z = z[:m] + 1j*z[m:]
    return (z[:p*n2].reshape(p,N,N), z[p*n2:2*p*n2].reshape(p,N,N),
            z[2*p*n2:(2*p+q)*n2].reshape(q,N,N), z[(2*p+q)*n2:].reshape(q,N,N))
def test(N, p, q, cplx, restarts, rng, tmax):
    T = crossing(N); n2 = N*N; m = (2*p+2*q)*n2; I = np.eye(N); best = 1e9; t0 = time.time(); r = 0
    def fun(z):
        F, Fp, G, Gp = unpack(z, N, p, q, cplx)
        d = (np.einsum('jxy,juv->xyuv', F, Fp) + np.einsum('kxu,kyv->xyuv', G, Gp) - T).ravel()
        return np.concatenate([d.real, d.imag]) if cplx else d
    def jac(z):
        F, Fp, G, Gp = unpack(z, N, p, q, cplx)
        J1 = np.einsum('xa,yb,juv->xyuvjab', I, I, Fp).reshape(n2*n2, -1)
        J2 = np.einsum('jxy,ua,vb->xyuvjab', F, I, I).reshape(n2*n2, -1)
        J3 = np.einsum('xa,ub,kyv->xyuvkab', I, I, Gp).reshape(n2*n2, -1)
        J4 = np.einsum('kxu,ya,vb->xyuvkab', G, I, I).reshape(n2*n2, -1)
        J = np.concatenate([J1, J2, J3, J4], axis=1)
        if not cplx: return J
        return np.block([[J.real, -J.imag], [J.imag, J.real]])
    while r < restarts and time.time()-t0 < tmax:
        r += 1
        z0 = rng.standard_normal(m*(2 if cplx else 1)) * rng.choice([0.3, 1.0])
        s = least_squares(fun, z0, jac=jac, method='trf', xtol=1e-13, ftol=1e-13, gtol=1e-13, max_nfev=250, x_scale='jac')
        best = min(best, np.sqrt(2*s.cost))
        if best < 1e-10: break
    return best, r
if __name__ == "__main__":
    rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
    tmax = float(sys.argv[3]) if len(sys.argv) > 3 else 25
    for spec in sys.argv[1].split(';'):
        N, p, q, f = spec.split(','); N, p, q = int(N), int(p), int(q)
        lbp = max(k*(N-k*q) for k in range(N+1)); lbq = max(k*(N-k*p) for k in range(N+1))
        best, r = test(N, p, q, f == 'c', 200, rng, tmax)
        print(f"N={N} (p,q)=({p},{q}) field={f.upper()} bound_ok={p>=lbp and q>=lbq} best_resid={best:.2e} restarts={r} ->", "FEASIBLE(num)" if best < 1e-10 else "not found", flush=True)
