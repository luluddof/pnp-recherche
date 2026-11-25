import sys, numpy as np
src = open('exact_points.py').read().split('for p in [1000003, 998244353]:')[0]
exec(src)
from math import factorial
kf = factorial(k); NN = N*N
Pn = {}
for m in shapes:
    P = np.zeros((N, N))
    for pm in perms: P[A[pm], np.arange(N)] += fdim[m]/kf*chi(m, pm)
    Pn[m] = P
K = np.zeros((NN, NN))
for m in shapes:
    for mp in shapes:
        if m != mp: K += np.kron(Pn[m], Pn[mp])
xs = np.arange(N); X, Y = np.meshgrid(xs, xs, indexing='ij'); X = X.ravel(); Y = Y.ravel(); col = X*N+Y
for m in shapes:
    M = np.zeros((NN, NN))
    for pi in perms:
        for sg in perms:
            M[A[pi][X]*N + A[sg][Y], col] += (fdim[m]/kf)**2*chi(m, compose(sg, pi))
    SM = np.zeros_like(M); SM[Y*N+X] = M; K += SM
S = np.zeros((NN, NN)); S[Y*N+X, X*N+Y] = 1
s = np.linalg.svd(S-K, compute_uv=False); r = (s > 1e-8).sum()
print("k", k, "float rank", r, "sv around gap", s[r-2:r+2])
C = K.reshape(N,N,N,N).transpose(0,2,1,3).reshape(NN,NN); s2 = np.linalg.svd(C, compute_uv=False); r2 = (s2 > 1e-8).sum()
print("osr", r2, s2[r2-1:r2+1], "K symmetric under pi(x)pi?", max(np.abs(K[np.ix_(A[pm][X]*N+A[pm][Y], A[pm][X]*N+A[pm][Y])]-K).max() for pm in perms))
