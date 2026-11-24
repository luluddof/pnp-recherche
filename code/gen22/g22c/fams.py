# Families b, c, d: generators of q signed (partial) permutation pairs (G_l, G'_l).
# Each family returns a list of (G, Gp) sparse int matrices and the convention C = sum G_l(x,u) G'_l(y,v).
import numpy as np, scipy.sparse as sp, itertools
from sympy import isprime, primitive_root

def P(N, perm, sign=None, dom=None):
    """M[perm[u], u] = sign[u] for u in dom (perm: array of length N)."""
    u = np.arange(N) if dom is None else np.asarray(dom)
    s = np.ones(len(u), dtype=np.int64) if sign is None else np.asarray(sign, dtype=np.int64)[u]
    return sp.csr_matrix((s, (np.asarray(perm)[u], u)), shape=(N, N))

# ---------- groups: elements 0..N-1, mult table ----------
def grp_cyclic(N):
    a = np.arange(N); return (a[:, None] + a[None, :]) % N, (-a) % N
def grp_elab(pr, n):
    N = pr ** n; el = np.array(list(itertools.product(range(pr), repeat=n)))
    idx = lambda v: int(np.ravel_multi_index(tuple(v % pr), (pr,) * n))
    T = np.array([[idx(el[i] + el[j]) for j in range(N)] for i in range(N)]); inv = np.array([idx(-el[i]) for i in range(N)])
    return T, inv
def grp_dihedral(m):   # order 2m: (r^i s^e), index i + m e
    N = 2 * m; T = np.zeros((N, N), int)
    for a in range(N):
        for b in range(N):
            i, e = a % m, a // m; j, f = b % m, b // m
            k = (i + (j if e == 0 else -j)) % m; T[a, b] = k + m * ((e + f) % 2)
    inv = np.array([int(np.where(T[a] == 0)[0][0]) for a in range(N)]); return T, inv
def grp_heis(pr):     # (x,y,z)(x',y',z') = (x+x', y+y', z+z'+x y')
    el = list(itertools.product(range(pr), repeat=3)); ix = {e: i for i, e in enumerate(el)}; N = len(el)
    T = np.array([[ix[((a[0] + b[0]) % pr, (a[1] + b[1]) % pr, (a[2] + b[2] + a[0] * b[1]) % pr)] for b in el] for a in el])
    inv = np.array([int(np.where(T[a] == 0)[0][0]) for a in range(N)]); return T, inv

def struct_slices(T, f, xs=None):
    """composition ansatz from a signed group table: G_x(a,c) = f(a,c)[a c = x], G'_x(b,d) = f(d,b)[d b = x]."""
    N = T.shape[0]; out = []
    for x in (range(N) if xs is None else xs):
        A, Cc = np.where(T == x)            # pairs (a,c) with a c = x
        G = sp.csr_matrix((f[A, Cc], (A, Cc)), shape=(N, N))
        Gp = sp.csr_matrix((f[A, Cc], (Cc, A)), shape=(N, N))   # G'(b,d) = f(d,b) [d b = x]: d = A, b = Cc
        out.append((G, Gp))
    return out

def translations(T, inv, gs, mode, sign=None):
    """G_g = left translation u -> g u (signed), G'_g = G_g (mode 'same'), G_g^T ('T'), or right translation ('R')."""
    N = T.shape[0]; out = []
    for g in gs:
        G = P(N, T[g, :], sign)
        if mode == 'same': Gp = G
        elif mode == 'T': Gp = G.T.tocsr()
        elif mode == 'R': Gp = P(N, T[:, g], sign)
        out.append((G, Gp))
    return out

# ---------- finite field F_r, r prime ----------
def chi(r):
    c = np.zeros(r, dtype=np.int64)
    sq = set((x * x) % r for x in range(1, r))
    for x in range(1, r): c[x] = 1 if x in sq else -1
    return c
def mult_family(r, A, mode, signmode):
    """G_a: x -> a x on F_r (0 fixed), signs: 'none', 'chi' (sign chi(u)), 'chia' (sign chi(a) const)."""
    c = chi(r); out = []
    for a in A:
        perm = (a * np.arange(r)) % r
        s = None if signmode == 'none' else (c if signmode == 'chi' else np.full(r, c[a]))
        G = P(r, perm, s)
        if mode == 'T': Gp = G.T.tocsr()
        elif mode == 'same': Gp = G
        elif mode == 'inv': Gp = P(r, (pow(int(a), r - 2, r) * np.arange(r)) % r, s)
        out.append((G, Gp))
    return out
def affine_family(r, pairs, mode, signmode):
    c = chi(r); out = []
    for a, b in pairs:
        perm = (a * np.arange(r) + b) % r
        s = None if signmode == 'none' else (c if signmode == 'chi' else np.full(r, c[a]))
        G = P(r, perm, s)
        Gp = G.T.tocsr() if mode == 'T' else G
        out.append((G, Gp))
    return out

# ---------- matching-shape families: G_l signed partial permutation, G'_l = eps G_l^T ----------
def matching_from_partition(N, parts, signs=None, eps=None):
    """parts: list of lists of pairs (x,u), each a partial permutation; G_l(x,u) = sign; G'_l = eps_l G_l^T."""
    out = []
    for l, M in enumerate(parts):
        if not M: continue
        xs = [m[0] for m in M]; us = [m[1] for m in M]
        s = [1] * len(M) if signs is None else signs[l]
        G = sp.csr_matrix((np.array(s, dtype=np.int64), (xs, us)), shape=(N, N))
        e = 1 if eps is None else eps[l]
        out.append((G, (e * G.T).tocsr()))
    return out
