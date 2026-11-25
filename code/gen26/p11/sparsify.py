# gauge-sparsify a numerical solution: G_l -> sum_k m_lk g G_k h^-1, G'_l -> sum_k (m^-T)_lk h G'_k g^-1  (C and X unchanged up to the
# symmetry x->g, v->g^-T, y->h, u->h^-T of X).  Minimise a smoothed L1 norm, then print.
import numpy as np, sys
from scipy.optimize import minimize
np.set_printoptions(precision=4, suppress=True, linewidth=200)
N, p, q = map(int, sys.argv[1].split(',')); z = np.load(sys.argv[4])
G = z[:q*N*N].reshape(q,N,N); Gp = z[q*N*N:].reshape(q,N,N)
def tr(w):
    g = w[:N*N].reshape(N,N); h = w[N*N:2*N*N].reshape(N,N); m = w[2*N*N:].reshape(q,q)
    gi, hi, mit = np.linalg.inv(g), np.linalg.inv(h), np.linalg.inv(m).T
    A = np.einsum('lk,kab->lab', m, np.einsum('xa,kab,bu->kxu', g, G, hi))
    B = np.einsum('lk,kab->lab', mit, np.einsum('ya,kab,bv->kyv', h, Gp, gi))
    return A, B
def obj(w, d):
    try: A, B = tr(w)
    except np.linalg.LinAlgError: return 1e9
    return np.sqrt(A**2 + d).sum() + np.sqrt(B**2 + d).sum()
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0); best = (1e9, None)
for r in range(int(sys.argv[3]) if len(sys.argv) > 3 else 12):
    w = np.concatenate([np.eye(N).ravel(), np.eye(N).ravel(), np.eye(q).ravel()]) + 0.5*rng.standard_normal(2*N*N+q*q)
    for d in (1e-2, 1e-4, 1e-6, 1e-8):
        w = minimize(obj, w, args=(d,), method='L-BFGS-B', options=dict(maxiter=3000)).x
    v = obj(w, 1e-10)
    if v < best[0]: best = (v, w.copy()); print("restart", r, "L1 =", v, flush=True)
A, B = tr(best[1]); np.save(sys.argv[4].replace('.npy','_sp.npy'), np.concatenate([A.ravel(), B.ravel()]))
print("G:"); print(A); print("G':"); print(B)
