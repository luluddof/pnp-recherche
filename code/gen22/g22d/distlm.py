# [N] distance^2 from T to J(p,q) over C by Levenberg-Marquardt on ALL factors: min ||T - sum_k a_k(xy) b_k(uv) - sum_l G_l(xu) H_l(yv)||^2
import numpy as np, sys, time
from scipy.optimize import least_squares
from dist import swap
def unpack(z, N, p, q):
    n = len(z)//2; w = z[:n] + 1j*z[n:]; s = N*N; o = 0
    a = w[o:o+p*s].reshape(p,N,N); o += p*s; b = w[o:o+p*s].reshape(p,N,N); o += p*s
    G = w[o:o+q*s].reshape(q,N,N); o += q*s; H = w[o:o+q*s].reshape(q,N,N)
    return a, b, G, H
def res(z, T, N, p, q):
    a, b, G, H = unpack(z, N, p, q)
    R = (T - np.einsum('kxy,kuv->xyuv', a, b) - np.einsum('kxu,kyv->xyuv', G, H)).ravel()
    return np.concatenate([R.real, R.imag])
def target(which, N, p, q, rng):
    if which == 'swap': return swap(N).astype(complex)
    if which == 'randJ':
        c = lambda sh: rng.standard_normal(sh) + 1j*rng.standard_normal(sh)
        T = np.einsum('kxy,kuv->xyuv', c((p,N,N)), c((p,N,N))) + np.einsum('kxu,kyv->xyuv', c((q,N,N)), c((q,N,N)))
        return T * N/np.linalg.norm(T)
if __name__ == "__main__":
    N, p, q, R, which = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    rng = np.random.default_rng(3); T = target(which, N, p, q, rng); out = []
    for r in range(R):
        z0 = rng.standard_normal(4*(p+q)*N*N) * 0.7
        s = least_squares(res, z0, args=(T,N,p,q), method='trf', max_nfev=4000, xtol=1e-15, ftol=1e-15, gtol=1e-15)
        v = 2*s.cost; out.append(v)
        print(f"  restart {r} dist^2={v:.4e} max|coef|={np.abs(s.x).max():.1f} nfev={s.nfev}", flush=True)
    print(f"N={N} (p,q)=({p},{q}) {which}: best {sorted(out)[:5]}", flush=True)
