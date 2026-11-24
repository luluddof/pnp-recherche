# G22c exact evaluator:  p = rank( d*F - sum_l kron(G_l, Gp_l) )   (Gp already multiplied by the common denominator d)
# F[(x,y),(u,v)] = [x=v][y=u];  kron(G,Gp)[(x,y),(u,v)] = G[x,u] Gp[y,v]  =>  C(x,y,u,v) = sum_l G_l(x,u) G'_l(y,v).
# Rank = sum over connected components of the bipartite row/col support graph; per block:
#   exact over Z (flint fmpz_mat) if block small, else mod two random ~62-bit primes (lower bounds on rank_Q; equal => reported).
import numpy as np, scipy.sparse as sp, flint, random
from scipy.sparse.csgraph import connected_components

PR = [4611686018427387847, 2305843009213693951, 2147483647]   # primes (2^62-57, 2^61-1, 2^31-1)

def build(N, terms, d=1):
    """terms: list of (G, Gp) integer arrays or scipy sparse. returns sparse int64 R (N^2 x N^2)."""
    idx = np.arange(N)
    rows = (idx[:, None] * N + idx[None, :]).ravel()       # (x,y)
    cols = (idx[None, :] * N + idx[:, None]).ravel()       # (u,v)=(y,x)
    R = sp.coo_matrix((np.full(N * N, d, dtype=np.int64), (rows, cols)), shape=(N * N, N * N)).tocsr()
    for G, Gp in terms:
        G = sp.csr_matrix(np.asarray(G, dtype=np.int64)) if not sp.issparse(G) else G.astype(np.int64)
        Gp = sp.csr_matrix(np.asarray(Gp, dtype=np.int64)) if not sp.issparse(Gp) else Gp.astype(np.int64)
        R = R - sp.kron(G, Gp, format='csr')
    R.eliminate_zeros()
    return R.tocsr()

def blocks(R):
    n, m = R.shape
    A = sp.bmat([[None, R != 0], [(R != 0).T, None]]).astype(np.int8)
    k, lab = connected_components(A, directed=False)
    out = []
    order = np.argsort(lab, kind='stable')
    bounds = np.searchsorted(lab[order], np.arange(k + 1))
    for c in range(k):
        nodes = order[bounds[c]:bounds[c + 1]]
        r = nodes[nodes < n]; cc = nodes[nodes >= n] - n
        if len(r) and len(cc):
            out.append((r, cc))
    return out

def rank_block(M, exact_limit=120, primes=PR[:2]):
    """M dense int64 ndarray. returns (rank, method)"""
    h, w = M.shape
    if min(h, w) == 0: return 0, 'e'
    if min(h, w) == 1: return int(np.any(M)), 'e'
    if max(h, w) <= exact_limit:
        return flint.fmpz_mat(M.tolist()).rank(), 'Q'
    rs = []
    for p in primes:
        rs.append(flint.nmod_mat((M % p).tolist(), p).rank())
    return max(rs), ('p' if len(set(rs)) == 1 else 'p!')

def rank(R, exact_limit=120, primes=PR[:2], verbose=False):
    tot = 0; meth = set(); sizes = []
    for r, c in blocks(R):
        B = R[r][:, c].toarray()
        rk, mth = rank_block(B, exact_limit, primes)
        tot += rk; meth.add(mth); sizes.append(len(r))
    if verbose: print('blocks', len(sizes), 'max', max(sizes) if sizes else 0, meth)
    return tot, ''.join(sorted(meth))

def p_of(N, terms, d=1, **kw):
    return rank(build(N, terms, d), **kw)

def dense_check(N, terms, d=1):
    """independent: full dense einsum build + flint fmpq rank over Q (whole matrix, no blocks)."""
    X = np.zeros((N, N, N, N), dtype=object)
    for x in range(N):
        for y in range(N): X[x, y, y, x] = d
    for G, Gp in terms:
        G = np.asarray(G.toarray() if sp.issparse(G) else G, dtype=object)
        Gp = np.asarray(Gp.toarray() if sp.issparse(Gp) else Gp, dtype=object)
        X = X - np.einsum('xu,yv->xyuv', G, Gp)
    M = X.reshape(N * N, N * N)
    return flint.fmpz_mat([[int(a) for a in row] for row in M]).rank()

def perm(N, f, s=None):
    """signed (partial) permutation: M[f(u)... ] convention M[x,u] = s(u) if x == f(u) (f(u) None => zero column)."""
    rows, cols, vals = [], [], []
    for u in range(N):
        x = f(u)
        if x is None: continue
        rows.append(x); cols.append(u); vals.append(1 if s is None else s(u))
    return sp.csr_matrix((np.array(vals, dtype=np.int64), (rows, cols)), shape=(N, N))

def build_M(N, terms):
    M = sp.csr_matrix((N * N, N * N), dtype=np.int64)
    for G, Gp in terms:
        G = sp.csr_matrix(np.asarray(G, dtype=np.int64)) if not sp.issparse(G) else G.astype(np.int64)
        Gp = sp.csr_matrix(np.asarray(Gp, dtype=np.int64)) if not sp.issparse(Gp) else Gp.astype(np.int64)
        M = M + sp.kron(G, Gp, format='csr')
    return M.tocsr()

def best_scale(N, terms, maxblock=2500, top=3):
    """min over scalars c of rank(F - c M), M = sum kron(G,Gp).  Candidates c = 1/lambda, lambda eigenvalue of F^T M
    (numeric, per block), ranked by numeric geometric multiplicity; the best rational candidates are re-evaluated exactly.
    returns list of (p_exact, lambda as (num,den))"""
    from fractions import Fraction
    M = build_M(N, terms)
    idx = np.arange(N)
    F = sp.coo_matrix((np.ones(N * N), ((idx[:, None] * N + idx[None, :]).ravel(), (idx[None, :] * N + idx[:, None]).ravel())), shape=(N * N, N * N)).tocsr()
    T = (F.T @ M).tocsr().astype(float)
    tot = {}
    S = (T != 0) + sp.eye(N * N, format='csr')
    k, lab = connected_components(S, directed=False)
    for c in range(k):
        nd = np.where(lab == c)[0]
        if len(nd) > maxblock: return None
        B = T[nd][:, nd].toarray()
        ev = np.linalg.eigvals(B)
        ev = ev[np.abs(ev.imag) < 1e-7].real
        cand = np.unique(np.round(ev, 6))
        for lam in cand:
            if abs(lam) < 1e-9: continue
            g = len(nd) - np.linalg.matrix_rank(B - lam * np.eye(len(nd)), tol=1e-7)
            if g > 0: tot[lam] = tot.get(lam, 0) + g
    best = sorted(tot.items(), key=lambda t: -t[1])[:top]
    out = []
    for lam, g in best:
        fr = Fraction(lam).limit_denominator(64)
        if abs(float(fr) - lam) > 1e-6: out.append((None, lam, N * N - g)); continue
        # rank(F - M/lam) = rank(num*F - den*M) with lam = num/den
        terms2 = [(G, (np.asarray(Gp.toarray() if sp.issparse(Gp) else Gp, dtype=np.int64) * fr.denominator)) for G, Gp in terms]
        r, mth = rank(build(N, terms2, d=fr.numerator))
        out.append((r, (fr.numerator, fr.denominator), N * N - g))
    return out
