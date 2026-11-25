# Sanity check of C3's site-pattern formulas: K = g(T) = sum_T g(T) prod_i (P_sym or P_anti)_i ; p = N^2 mu_b(g != parity),
# osr = sum_{U: c_U != 0} (n^2-1)^{|U|} in the centred basis tau_i = SWAP_i - 1/n.  Also Smolensky bound per instance.
import numpy as np, itertools
from math import comb
rng = np.random.default_rng(11)
def swap(n):
    S = np.zeros((n*n, n*n))
    for x in range(n):
        for y in range(n): S[y*n+x, x*n+y] = 1
    return S
def site_order(ops_per_site, n, k):
    # ops on (C^n x C^n)^{(x)k} site-major -> reorder to (C^n)^{(x)k} (x) (C^n)^{(x)k}
    M = ops_per_site[0]
    for o in ops_per_site[1:]: M = np.kron(M, o)
    dims = [n, n]*k
    T = M.reshape(dims + dims)
    perm = [2*i for i in range(k)] + [2*i+1 for i in range(k)]
    T = T.transpose(perm + [2*k + p for p in perm])
    return T.reshape(n**(2*k), n**(2*k))
for n, k in [(2,3), (2,4), (3,2), (3,3)]:
    N = n**k; Sw = swap(n); I = np.eye(n*n); Pa = (I - Sw)/2; Ps = (I + Sw)/2
    b = (n-1)/(2*n)
    for trial in range(4):
        # random g of degree <= D (multilinear in T), random coefficients, some agreeing with parity
        D = rng.integers(0, k+1)
        pts = list(itertools.product([0,1], repeat=k))
        # build via centred tau = 1 - 1/n - 2T ; choose random c_U for |U|<=D, some zero
        cU = {U: (rng.normal() if rng.random() < 0.7 else 0.0) for r in range(D+1) for U in itertools.combinations(range(k), r)}
        g = {T: sum(c*np.prod([1-1/n-2*T[i] for i in U]) for U, c in cU.items()) for T in pts}
        # force agreement with parity on a random subset by... (skip; generic g disagrees everywhere) -> use rounding trick
        K = sum(g[T]*site_order([Pa if T[i] else Ps for i in range(k)], n, k) for T in pts)
        S = site_order([Sw]*k, n, k)
        p = np.linalg.matrix_rank(S - K, tol=1e-8)
        mult = {T: np.prod([n*(n-1)/2 if t else n*(n+1)/2 for t in T]) for T in pts}
        pform = sum(mult[T] for T in pts if abs(g[T] - (-1)**sum(T)) > 1e-9)
        C = K.reshape(N,N,N,N).transpose(0,2,1,3).reshape(N*N,N*N)
        q = np.linalg.matrix_rank(C, tol=1e-8)
        qform = sum((n*n-1)**len(U) for U, c in cU.items() if c != 0)
        print(f"n={n} k={k} D={D}: p={p} formula {int(pform)}; osr={q} formula {qform}")
