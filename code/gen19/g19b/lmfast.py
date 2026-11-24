# Numerical feasibility of (N,p,q) over R: min || X - sum_j f_j(x,y) f'_j(u,v) - sum_l g_l(x,u) g'_l(y,v) ||.
# [N] only. Reports best residual and coefficient size (a residual -> 0 with exploding coefficients = border case).
import numpy as np, sys, time
from scipy.optimize import least_squares
def crossing(N):
    T = np.zeros((N,N,N,N))
    for x in range(N):
        for y in range(N): T[x,y,y,x] = 1.0
    return T
def unpack(z, N, p, q):
    n2 = N*N
    F = z[:p*n2].reshape(p,N,N); Fp = z[p*n2:2*p*n2].reshape(p,N,N)
    G = z[2*p*n2:(2*p+q)*n2].reshape(q,N,N); Gp = z[(2*p+q)*n2:].reshape(q,N,N)
    return F, Fp, G, Gp
def model(z, N, p, q):
    F, Fp, G, Gp = unpack(z, N, p, q)
    return np.einsum('jxy,juv->xyuv', F, Fp) + np.einsum('kxu,kyv->xyuv', G, Gp)
def jac(z, N, p, q):
    F, Fp, G, Gp = unpack(z, N, p, q); n2 = N*N; I = np.eye(N)
    J = np.zeros((N,N,N,N, (2*p+2*q)*n2)); o = 0
    # d/dF[j,a,b] = [x=a][y=b] Fp[j,u,v]
    JF = np.einsum('xa,yb,juv->xyuvjab', I, I, Fp).reshape(N,N,N,N,p*n2); J[..., o:o+p*n2] = JF; o += p*n2
    JFp = np.einsum('jxy,ua,vb->xyuvjab', F, I, I).reshape(N,N,N,N,p*n2); J[..., o:o+p*n2] = JFp; o += p*n2
    JG = np.einsum('xa,ub,kyv->xyuvkab', I, I, Gp).reshape(N,N,N,N,q*n2); J[..., o:o+q*n2] = JG; o += q*n2
    JGp = np.einsum('kxu,ya,vb->xyuvkab', G, I, I).reshape(N,N,N,N,q*n2); J[..., o:o+q*n2] = JGp
    return J.reshape(N**4, -1)
def solve(N, p, q, restarts=20, tmax=60, seed=0, tol=1e-11, verbose=False, nfev=400):
    T = crossing(N); rng = np.random.default_rng(seed); m = (2*p+2*q)*N*N; best = (1e9, None); t0 = time.time(); hist = []
    for r in range(restarts):
        z0 = rng.standard_normal(m) * rng.choice([0.2, 0.5, 1.0])
        try:
            s = least_squares(lambda z: (model(z,N,p,q)-T).ravel(), z0, jac=lambda z: jac(z,N,p,q), method='trf', x_scale='jac',
                          xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=nfev)
        except Exception as e:
            hist.append((9e9, -1)); continue
        res = np.sqrt(2*s.cost); hist.append((res, np.abs(s.x).max()))
        if res < best[0]: best = (res, s.x.copy())
        if verbose: print(f"   restart {r}: resid {res:.3e} max|z| {np.abs(s.x).max():.1f}", flush=True)
        if best[0] < tol or time.time() - t0 > tmax: break
    return best, hist
if __name__ == "__main__":
    for spec in sys.argv[1].split(';'):
        N, p, q = map(int, spec.split(','))
        tmax = float(sys.argv[2]) if len(sys.argv) > 2 else 60
        (res, z), hist = solve(N, p, q, restarts=200, tmax=tmax, seed=N*100+p*10+q, nfev=int(sys.argv[3]) if len(sys.argv) > 3 else 400)
        hs = sorted(hist)[:5]
        print(f"N={N} (p,q)=({p},{q}) best_resid={res:.3e} max|z|={np.abs(z).max():.2f} restarts={len(hist)} lowest5={[(float('%.2e'%a), float('%.0f'%b)) for a,b in hs]} ->",
              "FEASIBLE(num)" if res < 1e-10 else "not found", flush=True)
        if res < 1e-10: np.save(f"sol_{N}_{p}_{q}.npy", z)
