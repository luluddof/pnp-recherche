# Independent rebuild of C2's u=0 sector construction, d=2, N=2^k, exact mod primes.
# K = sum_{nu!=nu'} P_nu (x) P_nu'  +  sum_nu SWAP . M_nu,  M_nu = (f/k!)^2 sum_{pi,sig} chi_nu(sig pi) pi (x) sig
# (M_nu = SWAP_{S_nu} (x) 1_W on block (nu,nu)), P_nu = (f/k!) sum chi_nu(pi) pi.
import itertools, sys, numpy as np
from math import comb, factorial
from modrank import rank_mod
k = int(sys.argv[1]); d = 2; N = d**k
perms = list(itertools.permutations(range(k)))
idx = list(itertools.product(range(d), repeat=k))
pos = {t: i for i, t in enumerate(idx)}
def act(pm):  # vector v with (pi e)_x : e_x -> e_{v[x]}
    v = np.empty(N, dtype=np.int64)
    for i, t in enumerate(idx):
        new = [None]*k
        for s in range(k): new[pm[s]] = t[s]
        v[i] = pos[tuple(new)]
    return v
A = {pm: act(pm) for pm in perms}
def compose(a, b): return tuple(a[b[i]] for i in range(k))  # (a o b)
def cycles(pm):
    seen = set(); L = []
    for i in range(k):
        if i in seen: continue
        j = i; c = 0
        while j not in seen: seen.add(j); j = pm[j]; c += 1
        L.append(c)
    return L
def fixed_subsets(pm, m):
    dp = [1] + [0]*m
    for c in cycles(pm):
        for s in range(m, c-1, -1): dp[s] += dp[s-c]
    return dp[m]
def chi(m, pm):  # character of (k-m, m)
    return fixed_subsets(pm, m) - (fixed_subsets(pm, m-1) if m > 0 else 0)
shapes = list(range(k//2 + 1))
fdim = {m: comb(k, m) - (comb(k, m-1) if m else 0) for m in shapes}
wdim = {m: k - 2*m + 1 for m in shapes}
for p in [1000003, 998244353]:
    kf = factorial(k); inv = lambda x: pow(x % p, p-2, p)
    Pn = {}
    for m in shapes:
        P = np.zeros((N, N), dtype=np.int64)
        c0 = fdim[m] * inv(kf) % p
        for pm in perms:
            P[A[pm], np.arange(N)] = (P[A[pm], np.arange(N)] + c0 * (chi(m, pm) % p)) % p
        Pn[m] = P
    NN = N*N; K = np.zeros((NN, NN), dtype=np.int64)
    for m in shapes:
        for mp in shapes:
            if m != mp:
                K = (K + np.kron(Pn[m], Pn[mp]) % p) % p
    xs = np.arange(N); X, Y = np.meshgrid(xs, xs, indexing='ij'); X = X.ravel(); Y = Y.ravel()
    col = X*N + Y
    for m in shapes:
        M = np.zeros((NN, NN), dtype=np.int64)
        c0 = (fdim[m]**2) * inv(kf*kf) % p
        for pi in perms:
            for sg in perms:
                c = c0 * (chi(m, compose(sg, pi)) % p) % p
                if c == 0: continue
                row = A[pi][X]*N + A[sg][Y]
                M[row, col] = (M[row, col] + c) % p
        # SWAP . M : SWAP maps (x,y)->(y,x): row permutation
        sw = (Y*N + X)
        SM = np.zeros_like(M); SM[sw] = M[np.arange(NN)]  # (SWAP M)[(y,x),:] = M[(x,y),:]
        K = (K + SM) % p
    S = np.zeros((NN, NN), dtype=np.int64); S[Y*N+X, X*N+Y] = 1
    D = (S - K) % p
    prank = rank_mod(D, p)
    C = K.reshape(N, N, N, N).transpose(0, 2, 1, 3).reshape(NN, NN)
    q = rank_mod(C, p)
    pred = (N*N - sum(fdim[m]*wdim[m]**2 for m in shapes))//2
    print(f"k={k} N={N} prime={p}: rank(SWAP-K)={prank} (pred {pred}), osr={q}, dim commutant={comb(k+3,3)}", flush=True)
