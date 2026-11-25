# Does a (p,q)@N seed exist whose M_row = R S^-1 is NOT semisimple at 0 (Jordan block)?  N=2, (1,1): C = A(x)B, R = S - A(x)B
# rank 1; M = R S^-1 = I - (A(x)B) S, rank 1; nilpotent iff tr M = 0 iff tr(AB) = N^2 (tr((A(x)B)S) = tr(AB)).
# Numerical least squares over complex A,B: residual = singular values 2..4 of R, and tr(AB) - target.
import numpy as np, sys
from scipy.optimize import least_squares
N = int(sys.argv[1]); p = int(sys.argv[2]); q = int(sys.argv[3]); target = float(sys.argv[4]); trials = int(sys.argv[5])
S = np.zeros((N*N, N*N))
for x in range(N):
    for y in range(N): S[x*N+y, y*N+x] = 1
rng = np.random.default_rng(1)
def build(z):
    z = z[:len(z)//2] + 1j*z[len(z)//2:]
    A = z[:q*N*N].reshape(q, N, N); B = z[q*N*N:].reshape(q, N, N)
    K = sum(np.kron(A[l], B[l]) for l in range(q))       # operator (x,y)->(u,v): K[(x,y),(u,v)] = sum A[x,u]B[y,v]
    return K
def res(z):
    K = build(z); R = S - K; sv = np.linalg.svd(R, compute_uv=False)
    M = R @ S
    out = list(sv[p:]) + [abs(np.trace(M) - target)*1.0]
    return np.array(out)
best = None
for t in range(trials):
    z0 = rng.normal(size=2*2*q*N*N)
    r = least_squares(res, z0, xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=4000)
    if best is None or r.cost < best.cost: best = r
    if r.cost < 1e-20: break
K = build(best.x); R = S - K; M = R @ S
print("cost", best.cost, "sv(R)", np.round(np.linalg.svd(R, compute_uv=False), 8), "tr M", np.round(np.trace(M), 6))
print("eig M", np.round(np.linalg.eigvals(M), 6), " rank M^2:", np.linalg.matrix_rank(M @ M, 1e-6), " rank M:", np.linalg.matrix_rank(M, 1e-6))
