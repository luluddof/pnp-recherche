# Core tools for the full Bell model. Points z = c | d<<k in F2^{2k}; characters s = a | b<<k, chi_s(z) = (-1)^{|s&z|}.
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from fractions import Fraction
POP = None
def popc(x):
    x = np.asarray(x, dtype=np.int64); c = np.zeros_like(x)
    while np.any(x): c += x & 1; x = x >> 1
    return c
def beta(k):
    z = np.arange(4**k); c = z & ((1<<k)-1); d = z >> k
    return (1 - 2*(popc(c & d) & 1)).astype(np.int64)
def charmat(k, S):
    z = np.arange(4**k)
    S = np.asarray(S, dtype=np.int64)
    return (1 - 2*(popc(z[:,None] & S[None,:]) & 1)).astype(np.int8)
def fwht(a):
    a = a.astype(float).copy(); h = 1; n = len(a)
    while h < n:
        a = a.reshape(-1, 2*h); x = a[:, :h].copy(); y = a[:, h:].copy()
        a[:, :h] = x + y; a[:, h:] = x - y; a = a.reshape(n); h *= 2
    return a
def lad(M, b, w=None):
    """weighted L1 regression min sum w|Mx-b| via HiGHS. returns x."""
    n, q = M.shape
    if w is None: w = np.ones(n)
    c = np.concatenate([np.zeros(q), w, w])
    A = sp.hstack([sp.csr_matrix(M.astype(float)), sp.eye(n), -sp.eye(n)], format='csr')
    bounds = [(None, None)]*q + [(0, None)]*(2*n)
    r = linprog(c, A_eq=A, b_eq=b.astype(float), bounds=bounds, method='highs')
    return r.x[:q]
def rw_lad(M, b, iters=6, tau=1e-2):
    x = lad(M, b); best = (np.sum(np.abs(M@x - b) > 1e-7), x)
    for it in range(iters):
        res = np.abs(M@x - b); w = 1.0/(res + tau)
        x = lad(M, b, w); e = np.sum(np.abs(M@x - b) > 1e-7)
        if e < best[0]: best = (e, x)
    return best
def exact_errors(M, b, x, tol=1e-7):
    """Rigorous: take agreement set Z of float solution, pick q-rank rows, solve exactly over Q, count exact disagreements.
    returns (p_exact, nnz coefficients) for the exact solution (spectrum subset of S)."""
    n, q = M.shape
    Z = np.where(np.abs(M@x - b) <= tol)[0]
    rows = modp_rows(M[Z]); rows = Z[rows]
    # exact solve over Q on the column subspace: find columns basis too
    A = [[Fraction(int(v)) for v in M[r]] for r in rows]; rhs = [Fraction(int(b[r])) for r in rows]
    # Gaussian elimination rk x q, free vars set 0
    m = len(A); piv_cols = []; rI = 0
    A = [row[:] + [rhs[i]] for i, row in enumerate(A)]
    for col in range(q):
        pr = next((i for i in range(rI, m) if A[i][col] != 0), None)
        if pr is None: continue
        A[rI], A[pr] = A[pr], A[rI]; pv = A[rI][col]
        A[rI] = [v / pv for v in A[rI]]
        for i in range(m):
            if i != rI and A[i][col] != 0:
                f = A[i][col]; A[i] = [vi - f*vr for vi, vr in zip(A[i], A[rI])]
        piv_cols.append(col); rI += 1
        if rI == m: break
    xs = [Fraction(0)]*q
    for i, col in enumerate(piv_cols): xs[col] = A[i][q]
    D = 1
    for v in xs: D = D * v.denominator // np.gcd(D, v.denominator)
    X = np.array([int(v * D) for v in xs], dtype=object)
    Mo = M.astype(np.int64)
    if max(abs(int(v)) for v in X) * q < 2**62:
        val = Mo @ X.astype(np.int64)
    else:
        val = Mo.astype(object) @ X
    p = int(np.sum(val != b * D))
    return p, int(sum(1 for v in X if v != 0))

def modp_rows(A, P=2147483647):
    """indices of a maximal set of rows of integer matrix A independent mod P (hence over Q)."""
    A = (A.astype(np.int64) % P).copy(); n, q = A.shape
    alive = np.ones(n, bool); chosen = []
    for col in range(q):
        cand = np.where(alive & (A[:, col] != 0))[0]
        if len(cand) == 0: continue
        r = cand[0]; chosen.append(r); alive[r] = False
        inv = pow(int(A[r, col]), P-2, P)
        piv = (A[r] * inv) % P
        f = A[:, col].copy(); f[r] = 0
        idx = np.where(f != 0)[0]
        for i0 in range(0, len(idx), 4096):
            ii = idx[i0:i0+4096]
            A[ii] = (A[ii] - (f[ii, None] * piv[None, :]) % P) % P
    return np.array(chosen, dtype=np.int64)
