# shared helpers: random points of J(p,q), of Z3 = S_(14|23), the swap; all mod a prime.
import numpy as np, random
from flint import nmod_mat
PR = 32003
def rnd(shape, rng): return rng.integers(0, PR, size=shape, dtype=np.int64)
def swap(N):
    X = np.zeros((N,)*4, dtype=np.int64)
    for x in range(N):
        for y in range(N): X[x,y,y,x] = 1
    return X
def ptJ(N, p, q, rng):
    """A + B: rank_(xy|uv) A <= p, rank_(xu|yv) B <= q; returns X[x,y,u,v] mod PR"""
    a = rnd((p,N,N), rng); b = rnd((p,N,N), rng); G = rnd((q,N,N), rng); H = rnd((q,N,N), rng)
    A = np.einsum('kxy,kuv->xyuv', a, b) % PR
    B = np.einsum('kxu,kyv->xyuv', G, H) % PR
    return (A + B) % PR
def ptZ3(N, rng):
    M = rnd((N,N), rng); Mp = rnd((N,N), rng)
    return np.einsum('xv,yu->xyuv', M, Mp) % PR
def rank(M):
    M = np.asarray(M, dtype=np.int64) % PR
    if M.size == 0: return 0
    return nmod_mat(M.shape[0], M.shape[1], [int(v) for v in M.ravel()], PR).rank()
