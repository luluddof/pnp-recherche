# [N] adversarial test of the refined p=1 bounds (Theorem A'):
#  kappa(v,B) <= (1+y)^2/(1-w1)   and   kappa(v,B) <= 1 + 2(1+y)^2/(1-w2),
#  y = sigma1*sigma2 (Schmidt coeffs of v), w1 >= w2 the two largest |<v|z_j>|^2 over the product right-singular basis of B.
import numpy as np, sys
from scipy.optimize import minimize
N, R = int(sys.argv[1]), int(sys.argv[2]); n = N*N
Fm = np.zeros((n, n))
for x in range(N):
    for y in range(N): Fm[x*N+y, y*N+x] = 1
def parts(z):
    w = z[:len(z)//2] + 1j*z[len(z)//2:]
    G = w[:n].reshape(N, N); Gp = w[n:2*n].reshape(N, N); v = w[2*n:]; v = v/np.linalg.norm(v)
    return G, Gp, v
def ratio(z):
    G, Gp, v = parts(z); B = np.kron(G, Gp); Q = np.eye(n) - np.outer(v, v.conj())
    kap = abs(np.trace(Fm @ B @ Q))**2/np.linalg.norm(B @ Q)**2
    sg = np.linalg.svd(v.reshape(N, N), compute_uv=False); y = sg[0]*sg[1]
    _, _, v1 = np.linalg.svd(G); _, _, v2 = np.linalg.svd(Gp)
    ov = sorted([abs(np.vdot(np.kron(v1[i].conj(), v2[k].conj()), v))**2 for i in range(N) for k in range(N)], reverse=True)
    b1 = (1+y)**2/(1-ov[0]); b2 = 1 + 2*(1+y)**2/(1-ov[1])
    return kap/min(b1, b2), kap
rng = np.random.default_rng(3); worst = 0; kmax = 0
for r in range(R):
    s = minimize(lambda z: -ratio(z)[0], rng.standard_normal(6*n), method='Nelder-Mead', options=dict(maxiter=20000, maxfev=20000, xatol=1e-10, fatol=1e-12))
    rr, kk = ratio(s.x); worst = max(worst, rr); kmax = max(kmax, kk)
    print(r, "kappa/bound = %.6f  kappa = %.5f" % (rr, kk), flush=True)
print("RESULT N=%d worst kappa/bound = %.6f (must be <= 1)" % (N, worst))
