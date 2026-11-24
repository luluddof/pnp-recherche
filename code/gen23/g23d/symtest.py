# [N] test of the "symmetric compression" strengthening: sup over Y of osr <= q of |<Omega|Y|Omega>|^2 / ||P_S Y P_S||^2
import numpy as np, sys
from scipy.optimize import minimize
N, q, R = map(int, sys.argv[1:4])
Om = np.eye(N).reshape(N*N)
Fm = np.zeros((N*N, N*N))
for x in range(N):
    for y in range(N): Fm[x*N+y, y*N+x] = 1
PS = (np.eye(N*N) + Fm)/2
rng = np.random.default_rng(0); best = 0
for r in range(R):
    def f(z):
        w = z[:len(z)//2] + 1j*z[len(z)//2:]; w = w.reshape(q, 2, N, N)
        Y = sum(np.kron(w[l,0], w[l,1]) for l in range(q))
        num = abs(Om @ Y @ Om)**2; den = np.linalg.norm(PS @ Y @ PS)**2
        return -num/den
    s = minimize(f, rng.standard_normal(4*q*N*N), method='BFGS', options=dict(gtol=1e-10, maxiter=5000))
    best = max(best, -s.fun); print(r, -s.fun, flush=True)
print("RESULT", N, q, best)
