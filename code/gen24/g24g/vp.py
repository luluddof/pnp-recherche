# [N] Variable-projection solver for d^2(N;p,q) = inf ||F - A - B||^2, rank A <= p (Fl1 = xy|uv), osr B <= q (Fl2 = xu|yv).
# Eliminates A exactly (Eckart-Young on F - B in Fl1); optimises B = sum_l G_l (x) H_l (complex) with L-BFGS.
# f(B) = sum_{i>p} sigma_i(F-B)^2, gradient by envelope theorem: d f = -2 Re <E, dB>, E = residual.
# By symmetry c(p,q) = c(q,p) [realignment], we always optimise over the side with the smaller parameter... (here: B side = q).
import numpy as np, sys, time, json
from scipy.optimize import minimize

def swap(N):
    X = np.zeros((N,)*4)
    for x in range(N):
        for y in range(N): X[x,y,y,x] = 1
    return X

def fg(z, F, N, p, q):
    n = len(z)//2; w = z[:n] + 1j*z[n:]; s = N*N
    G = w[:q*s].reshape(q,N,N); H = w[q*s:].reshape(q,N,N)
    B = np.einsum('kxu,kyv->xyuv', G, H)
    M = (F - B).reshape(s, s)                      # Fl1
    U, S, Vh = np.linalg.svd(M)
    E = (U[:, p:] * S[p:]) @ Vh[p:, :]
    f = float(np.sum(S[p:]**2))
    E = E.reshape(N,N,N,N)
    gG = -2*np.einsum('xyuv,kyv->kxu', E, H.conj())
    gH = -2*np.einsum('xyuv,kxu->kyv', E, G.conj())
    g = np.concatenate([gG.ravel(), gH.ravel()])
    return f, np.concatenate([g.real, g.imag]), S

def run(N, p, q, R, seed=1, maxiter=20000, scale=None, verbose=True, init=None):
    if q > p: p, q = q, p                          # optimise over the smaller side (symmetry c(p,q)=c(q,p))
    F = swap(N).astype(complex); rng = np.random.default_rng(seed); out = []
    for r in range(R):
        sc = scale if scale else rng.choice([0.3, 0.7, 1.0, 1.5])
        z0 = rng.standard_normal(4*q*N*N) * sc
        if init is not None: z0 = init + (r > 0) * 1e-3 * rng.standard_normal(len(init)) * np.abs(init).mean()
        res = minimize(lambda z: fg(z, F, N, p, q)[:2], z0, jac=True, method='L-BFGS-B',
                       options=dict(maxiter=maxiter, maxfun=2*maxiter, ftol=1e-16, gtol=1e-12, maxcor=30))
        f, g, S = fg(res.x, F, N, p, q)
        if not out or f < min(o['d2'] for o in out): np.save(f"best_{N}_{p}_{q}_s{seed}.npy", res.x)
        out.append(dict(d2=f, coef=float(np.abs(res.x).max()), nit=res.nit, sp=[round(float(t),4) for t in S[max(0,p-2):p+4]]))
        if verbose: print(f" r{r} d2={f:.6f} c={N*N-f:.4f} coef={out[-1]['coef']:.1f} nit={res.nit} sv[p-2:p+4]={out[-1]['sp']}", flush=True)
    ds = sorted(o['d2'] for o in out)
    return ds, out

if __name__ == "__main__":
    N, p, q, R = map(int, sys.argv[1:5]); seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    maxiter = int(sys.argv[6]) if len(sys.argv) > 6 else 20000
    init = np.load(sys.argv[7]) if len(sys.argv) > 7 else None
    t0 = time.time(); ds, out = run(N, p, q, R, seed, maxiter, init=init)
    print(f"RESULT N={N} p={p} q={q} best_d2={ds[0]:.6f} c={N*N-ds[0]:.4f} sorted_d2={[round(d,4) for d in ds]} t={time.time()-t0:.0f}s", flush=True)
    json.dump(out, open(f"vp_{N}_{p}_{q}_s{seed}.json", "w"))
