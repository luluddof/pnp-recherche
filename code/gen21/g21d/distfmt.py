# [N] LM distance^2 to J(p,q) in a rectangular format (nx,ny,nu,nv): A = sum_k a_k(x,y) b_k(u,v), B = sum_l G_l(x,u) H_l(y,v)
import numpy as np, sys
from scipy.optimize import least_squares
def build(z, f, p, q):
    nx, ny, nu, nv = f; n = len(z)//2; w = z[:n] + 1j*z[n:]; o = 0; out = []
    for cnt, sh in [(p,(nx,ny)),(p,(nu,nv)),(q,(nx,nu)),(q,(ny,nv))]:
        s = cnt*sh[0]*sh[1]; out.append(w[o:o+s].reshape(cnt, *sh)); o += s
    a, b, G, H = out
    return np.einsum('kxy,kuv->xyuv', a, b) + np.einsum('kxu,kyv->xyuv', G, H)
def nparam(f, p, q):
    nx, ny, nu, nv = f; return 2*(p*(nx*ny+nu*nv) + q*(nx*nu+ny*nv))
def res(z, T, f, p, q):
    R = (T - build(z, f, p, q)).ravel(); return np.concatenate([R.real, R.imag])
if __name__ == "__main__":
    f = tuple(map(int, sys.argv[1].split(','))); p, q, R, which = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    rng = np.random.default_rng(5)
    if which == 'swapproj':     # delta_{xv} delta_{yu} restricted to the format
        T = np.zeros(f, dtype=complex)
        for x in range(f[0]):
            for y in range(f[1]):
                if x < f[3] and y < f[2]: T[x,y,y,x] = 1
    else:
        T = build(rng.standard_normal(nparam(f,p,q)), f, p, q); T *= np.sqrt(6)/np.linalg.norm(T)
    out = []
    for r in range(R):
        s = least_squares(res, rng.standard_normal(nparam(f,p,q))*0.7, args=(T,f,p,q), method='trf', max_nfev=3000, xtol=1e-15, ftol=1e-15, gtol=1e-15)
        out.append(2*s.cost); print(f"  restart {r} dist^2={2*s.cost:.4e} max|coef|={np.abs(s.x).max():.1f}", flush=True)
    print(f"format {f} (p,q)=({p},{q}) {which} ||T||^2={np.linalg.norm(T)**2:.1f}: best {['%.3e'%v for v in sorted(out)[:6]]}", flush=True)
