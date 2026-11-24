# Rank-method evaluator: Phi(X) = sum_t A_t Fl_{s_t}(X) B_t, exact ranks mod a large prime (python-flint nmod_mat).
# c1 = generic rank of Phi on rank-one operators, c2 = on osr-one operators (generic = max, by semicontinuity;
# computed at 2 random points, max taken), R = rank Phi(F).
import numpy as np, sys, random
from flint import nmod_mat
PR = 1000003
def rank(M):
    M = np.asarray(M) % PR
    return nmod_mat(M.shape[0], M.shape[1], [int(v) for v in M.flatten()], PR).rank()
def fl(X4, s, N):
    # X4 indices (x,y,u,v)
    if s == 1: return X4.reshape(N*N, N*N)                       # (xy|uv)
    if s == 2: return X4.transpose(0,2,1,3).reshape(N*N, N*N)    # (xu|yv)
    if s == 3: return X4.transpose(0,3,1,2).reshape(N*N, N*N)    # (xv|yu)
def swap(N):
    X = np.zeros((N,N,N,N), dtype=np.int64)
    for x in range(N):
        for y in range(N): X[x,y,y,x] = 1
    return X
def apply(phi, X4, N):
    out = None
    for (s, A, B) in phi:
        T = (A @ (fl(X4, s, N) % PR) % PR) @ B % PR
        out = T if out is None else (out + T) % PR
    return out
def rnd(*sh): return np.random.randint(0, PR, size=sh).astype(object)
def rank1op(N):
    u = rnd(N*N); v = rnd(N*N); return np.outer(u, v).reshape(N,N,N,N) % PR
def osr1op(N):
    G = rnd(N,N); H = rnd(N,N)   # X(x,y,u,v) = G(x,u) H(y,v)
    return np.einsum('xu,yv->xyuv', G, H) % PR
def profile(phi, N, reps=2):
    c1 = max(rank(apply(phi, rank1op(N), N)) for _ in range(reps))
    c2 = max(rank(apply(phi, osr1op(N), N)) for _ in range(reps))
    R = rank(apply(phi, swap(N).astype(object), N))
    return c1, c2, R
