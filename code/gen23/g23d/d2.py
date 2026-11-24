# [N] d^2(N;p,q) = inf ||F - A - B||_F^2, rank A <= p (Fl1 = xy|uv), osr B <= q (Fl2 = xu|yv).
# Levenberg-Marquardt on all factors (complex), analytic Jacobian.  Prints residual singular values in Fl1/Fl2/Fl3.
import numpy as np, sys, json, time
from scipy.optimize import least_squares

def swap(N):
    X = np.zeros((N,)*4, complex)
    for x in range(N):
        for y in range(N): X[x,y,y,x] = 1
    return X

def unpack(z, N, p, q):
    n = len(z)//2; w = z[:n] + 1j*z[n:]; s = N*N; o = 0
    a = w[o:o+p*s].reshape(p,N,N); o += p*s; b = w[o:o+p*s].reshape(p,N,N); o += p*s
    G = w[o:o+q*s].reshape(q,N,N); o += q*s; H = w[o:o+q*s].reshape(q,N,N)
    return a, b, G, H

def model(z, N, p, q):
    a, b, G, H = unpack(z, N, p, q)
    return np.einsum('kxy,kuv->xyuv', a, b), np.einsum('kxu,kyv->xyuv', G, H)

def res(z, T, N, p, q):
    A, B = model(z, N, p, q); R = (T - A - B).ravel()
    return np.concatenate([R.real, R.imag])

def jac(z, T, N, p, q):
    a, b, G, H = unpack(z, N, p, q); I = np.eye(N)
    Ja = -np.einsum('xX,yY,kuv->xyuvkXY', I, I, b).reshape(N**4, -1)
    Jb = -np.einsum('kxy,uU,vV->xyuvkUV', a, I, I).reshape(N**4, -1)
    JG = -np.einsum('xX,uU,kyv->xyuvkXU', I, I, H).reshape(N**4, -1)
    JH = -np.einsum('kxu,yY,vV->xyuvkYV', G, I, I).reshape(N**4, -1)
    J = np.concatenate([Ja, Jb, JG, JH], axis=1)
    return np.block([[J.real, -J.imag], [J.imag, J.real]])

def sv(E, perm, N):
    return np.linalg.svd(E.transpose(perm).reshape(N*N, N*N), compute_uv=False)

def solve(T, N, p, q, z0, nfev=6000):
    s = least_squares(res, z0, jac=jac, args=(T,N,p,q), method='lm', max_nfev=nfev, xtol=1e-15, ftol=1e-15, gtol=1e-15)
    A, B = model(s.x, N, p, q); E = T - A - B
    return float(np.linalg.norm(E)**2), s.x, E, A, B

def target(which, N, p, q, rng):
    if which == 'swap': return swap(N)
    c = lambda sh: rng.standard_normal(sh) + 1j*rng.standard_normal(sh)
    T = np.einsum('kxy,kuv->xyuv', c((p,N,N)), c((p,N,N))) + np.einsum('kxu,kyv->xyuv', c((q,N,N)), c((q,N,N)))
    return T * N/np.linalg.norm(T)

if __name__ == "__main__":
    N, p, q, R = map(int, sys.argv[1:5]); which = sys.argv[5] if len(sys.argv) > 5 else 'swap'
    nfev = int(sys.argv[6]) if len(sys.argv) > 6 else 6000
    seed = int(sys.argv[7]) if len(sys.argv) > 7 else 1
    rng = np.random.default_rng(seed); T = target(which, N, p, q, rng); out = []; t0 = time.time()
    for r in range(R):
        z0 = rng.standard_normal(4*(p+q)*N*N) * rng.choice([0.3, 0.7, 1.5])
        v, x, E, A, B = solve(T, N, p, q, z0, nfev)
        s1, s2, s3 = sv(E,(0,1,2,3),N), sv(E,(0,2,1,3),N), sv(E,(0,3,1,2),N)
        f = lambda s: [round(float(t),4) for t in s if t > 1e-3]
        out.append(dict(d2=v, coef=float(np.abs(x).max()), fl1=f(s1), fl2=f(s2), fl3=f(s3), opn=float(s1[0])))
        print(f" r{r} d2={v:.6f} coef={out[-1]['coef']:.1f} Fl1sv={f(s1)} Fl2sv={f(s2)} Fl3sv={f(s3)}", flush=True)
    ds = sorted(o['d2'] for o in out)
    print(f"RESULT N={N} p={p} q={q} {which} best={ds[0]:.6f} sorted={[round(d,4) for d in ds]} t={time.time()-t0:.0f}s", flush=True)
    json.dump(out, open(f"out_{which}_{N}_{p}_{q}.json", "w"))
