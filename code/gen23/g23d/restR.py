# [N] conjecture R: B_r = G(x)G' minus its top singular term; sup_v |tr(F B_r Q)|^2/||B_r Q||^2 (Q = 1 - vv^*) <= 2 ?
import numpy as np, sys
from scipy.optimize import minimize
N, R = int(sys.argv[1]), int(sys.argv[2])
Fm = np.zeros((N*N, N*N))
for x in range(N):
    for y in range(N): Fm[x*N+y, y*N+x] = 1
def val(z):
    n = N*N; w = z[:len(z)//2] + 1j*z[len(z)//2:]
    G = w[:n].reshape(N, N); Gp = w[n:2*n].reshape(N, N); v = w[2*n:]; v = v/np.linalg.norm(v)
    B = np.kron(G, Gp)
    U, s, Vh = np.linalg.svd(B); Br = B - s[0]*np.outer(U[:, 0], Vh[0])
    Q = np.eye(n) - np.outer(v, v.conj())
    BQ = Br @ Q
    return abs(np.trace(Fm @ BQ))**2/np.linalg.norm(BQ)**2
rng = np.random.default_rng(1); best = 0
for r in range(R):
    s = minimize(lambda z: -val(z), rng.standard_normal(2*(3*N*N)), method='Nelder-Mead' if False else 'BFGS', options=dict(maxiter=4000))
    best = max(best, -s.fun); print(r, -s.fun, flush=True)
print("RESULT", N, best)
