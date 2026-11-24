# [N] distance^2 from a target T to J(p,q) over C:  min over complex B = sum_{l<=q} G_l(x)H_l (Fl2) of
#   tail_p(T - B) = sum_{i>p} sigma_i(Fl1(T - B))^2.   L-BFGS on real+imag parts, random restarts.
import numpy as np, sys, time
from scipy.optimize import minimize
def swap(N):
    X = np.zeros((N,)*4)
    for x in range(N):
        for y in range(N): X[x,y,y,x] = 1
    return X
def fg(z, N, p, q, T):
    n = q*N*N
    w = z[:2*n] + 1j*z[2*n:]
    G = w[:n].reshape(q,N,N); H = w[n:].reshape(q,N,N)
    B = np.einsum('kxu,kyv->xyuv', G, H)
    M = (T - B).reshape(N*N, N*N)
    U, s, Vt = np.linalg.svd(M)
    Et = ((U[:, p:] * s[p:]) @ Vt[p:]).reshape(N,N,N,N)
    gG = np.einsum('xyuv,kyv->kxu', Et.conj(), H); gH = np.einsum('xyuv,kxu->kyv', Et.conj(), G)
    g = np.concatenate([gG.ravel(), gH.ravel()])
    return float((s[p:]**2).sum()), np.concatenate([-2*g.real, 2*g.imag])
def run(T, N, p, q, restarts, seed, maxit=int(__import__("os").environ.get("MAXIT","3000"))):
    rng = np.random.default_rng(seed); out = []
    for r in range(restarts):
        z0 = rng.standard_normal(4*q*N*N) * rng.choice([0.3, 1.0, 2.0])
        res = minimize(fg, z0, args=(N,p,q,T), jac=True, method='L-BFGS-B', options=dict(maxiter=maxit, maxfun=2*maxit, ftol=1e-30, gtol=1e-13))
        out.append((res.fun, float(np.abs(res.x).max())))
        print('   restart', r, '%.4e' % res.fun, 'max|coef| %.1f' % out[-1][1], 'nit', res.nit, flush=True)
    return sorted(out)
if __name__ == "__main__":
    N, p, q, R = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    which = sys.argv[5]
    rng = np.random.default_rng(7)
    if which == 'swap': T = swap(N)
    elif which == 'randJ':   # random complex point of J(p,q) scaled to ||T||^2 = N^2 (calibration: must give 0)
        a = rng.standard_normal((p,N,N))+1j*rng.standard_normal((p,N,N)); b = rng.standard_normal((p,N,N))+1j*rng.standard_normal((p,N,N))
        G = rng.standard_normal((q,N,N))+1j*rng.standard_normal((q,N,N)); H = rng.standard_normal((q,N,N))+1j*rng.standard_normal((q,N,N))
        T = np.einsum('kxy,kuv->xyuv', a, b) + np.einsum('kxu,kyv->xyuv', G, H); T *= N/np.linalg.norm(T)
    elif which == 'generic':  # random operator, norm N (calibration: positive)
        T = rng.standard_normal((N,)*4) + 1j*rng.standard_normal((N,)*4); T *= N/np.linalg.norm(T)
    t0 = time.time(); out = run(T, N, p, q, R, 11)
    print(f"N={N} (p,q)=({p},{q}) target={which} ||T||^2={np.linalg.norm(T)**2:.1f} restarts={R} best tails (value, max|coef|): {[(float('%.3e'%a), round(b,1)) for a,b in out[:8]]} time={time.time()-t0:.0f}s", flush=True)
