"""Numerical check of Appendix B pieces (Heisenberg/Fock construction).
Pieces A_s(X,Y,U,V) = [t^U] Q_s(a_X, b_V) t^Y over partitions X,Y,U,V of h.
Checks: (1) orig-order and exchanged-order Q_s agree (Lemma 9a),
        (2) sum_s A_s = h! [X=V][Y=U] (Lemma 9c + Lemma 10),
        (3) exact row/column ranks of each A_s, and resulting bound on P(p(h)).
"""
import sys
from functools import lru_cache
from math import comb, factorial
from collections import defaultdict
import flint

def partitions(n, maxp=None):
    if maxp is None: maxp = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxp), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest

@lru_cache(None)
def parts_list(E):
    return list(partitions(E)) if E >= 0 else []

@lru_cache(None)
def p(E):
    return len(parts_list(E))

@lru_cache(None)
def index(E):
    return {mu: i for i, mu in enumerate(parts_list(E))}

def gbinom(k, j):
    if j < 0: return 0
    num = 1
    for t in range(j): num *= (k - t)
    return num // factorial(j)

def add_part(sigma, l):
    return tuple(sorted(sigma + (l,), reverse=True))

def rem_part(sigma, k):
    s = list(sigma); s.remove(k); return tuple(s)

@lru_cache(None)
def apply_mono(parts, r, sigma):
    """(t^parts)[r] applied to t^sigma, via recursion (O2). Returns tuple of (mono, coef)."""
    if not parts:
        return ((sigma, 1),) if r == 0 else ()
    E = sum(sigma)
    if E + r < 0: return ()
    n, rest = parts[0], parts[1:]
    out = defaultdict(int)
    for l in range(n, E + r + 1):
        c = comb(l - 1, n - 1)
        for m, v in apply_mono(rest, r - l, sigma):
            out[add_part(m, l)] += c * v
    for k in set(sigma):
        c = gbinom(-k - 1, n - 1) * k * sigma.count(k)
        for m, v in apply_mono(rest, r + k, rem_part(sigma, k)):
            out[m] += c * v
    return tuple((m, v) for m, v in out.items() if v)

# direct implementation of (47) for cross-check
def apply_direct(parts, r, sigma):
    E = sum(sigma)
    out = defaultdict(int)
    m = len(parts)
    rng = [k for k in range(-(E + abs(r) + 5), E + 1) if k != 0]
    def rec(i, ks):
        if i == m:
            if sum(ks) != -r: return
            c = 1
            for mu, k in zip(parts, ks): c *= gbinom(-k - 1, mu - 1)
            if c == 0: return
            f = {sigma: 1}
            for k in ks:
                if k > 0:
                    g = defaultdict(int)
                    for mono, v in f.items():
                        if k in mono: g[rem_part(mono, k)] += v * k * mono.count(k)
                    f = g
            for k in ks:
                if k < 0:
                    f = {add_part(mono, -k): v for mono, v in f.items()}
            for mono, v in f.items(): out[mono] += c * v
            return
        for k in rng: rec(i + 1, ks + [k])
    if m == 0:
        return {sigma: 1} if r == 0 else {}
    rec(0, [])
    return {k: v for k, v in out.items() if v}

def poly_mul(f, g):
    out = defaultdict(int)
    for a, x in f.items():
        for b, y in g.items():
            out[tuple(sorted(a + b, reverse=True))] += x * y
    return {k: v for k, v in out.items() if v}

def op_matrix(poly, r, E):
    """Matrix of poly[r]: V_E -> V_{E+r}, as list of rows (dict col->val)."""
    rows, cols = p(E + r), p(E)
    M = [[0] * cols for _ in range(rows)]
    if rows == 0 or cols == 0: return M
    idx = index(E + r)
    for j, sigma in enumerate(parts_list(E)):
        for mono, c in poly.items():
            for m, v in apply_mono(mono, r, sigma):
                M[idx[m]][j] += c * v
    return M

def matmul(A, B):
    n, k = len(A), len(B)
    m = len(B[0]) if k else 0
    if k == 0: return [[0] * m for _ in range(n)]
    return [[sum(A[i][t] * B[t][j] for t in range(k)) for j in range(m)] for i in range(n)]

def run(h):
    P = parts_list(h); M = len(P)
    y = {}
    for n in range(1, h + 1):
        y[n] = {(j,): comb(n, j) for j in range(1, n + 1)}
    def zee(mu):
        z = 1
        for k in set(mu):
            mk = mu.count(k); z *= k ** mk * factorial(mk)
        return z
    a = {}
    for mu in P:
        f = {(): 1}
        for part in mu:
            f = poly_mul(f, {k: (-1) ** (part - 1) * v for k, v in y[part].items()})
        a[mu] = f
    b = {nu: {nu: factorial(h) // zee(nu)} for nu in P}
    # a<i> = sum_{h'} C(h'+h, i) abar_{h'}
    def a_i(mu, i):
        return {k: v * comb(sum(k) + h, i) for k, v in a[mu].items() if v * comb(sum(k) + h, i)}
    S = range(-2 * h, 2 * h + 1)
    A = {}
    lemma9a_ok = True
    for s in S:
        T = [[[[0] * M for _ in range(M)] for _ in range(M)] for _ in range(M)]  # T[X][Y][U][V]
        for X in range(M):
            for V in range(M):
                orig = [[0] * M for _ in range(M)]  # [U][Y]
                exch = [[0] * M for _ in range(M)]
                for i in range(2 * h + 1):
                    ai = a_i(P[X], i); sg = (-1) ** i
                    # orig: b[h-s-i]: V_h -> V_{2h-s-i}, then a<i>[s-h+i]: V_{2h-s-i} -> V_h
                    w = 2 * h - s - i
                    if w >= 0:
                        Mb = op_matrix(b[P[V]], h - s - i, h)
                        Ma = op_matrix(ai, s - h + i, w)
                        prod = matmul(Ma, Mb)
                        for U in range(M):
                            for Y in range(M): orig[U][Y] += sg * prod[U][Y]
                    w2 = s + i
                    if w2 >= 0:
                        Ma = op_matrix(ai, s - h + i, h)
                        Mb = op_matrix(b[P[V]], h - s - i, w2)
                        prod = matmul(Mb, Ma)
                        for U in range(M):
                            for Y in range(M): exch[U][Y] += sg * prod[U][Y]
                if orig != exch: lemma9a_ok = False
                for U in range(M):
                    for Y in range(M): T[X][Y][U][V] = orig[U][Y]
        A[s] = T
    # identity check
    hf = factorial(h)
    ident_ok = True
    for X in range(M):
        for Y in range(M):
            for U in range(M):
                for V in range(M):
                    tot = sum(A[s][X][Y][U][V] for s in S)
                    if tot != hf * (X == V) * (Y == U): ident_ok = False
    print(f"h={h} p(h)={M} lemma9a(both orders equal)={lemma9a_ok} sum_s A_s == h!*crossing: {ident_ok}")
    total_min = 0
    rows_out = []
    for s in S:
        T = A[s]
        Rrow = flint.fmpz_mat([[T[X][Y][U][V] for U in range(M) for V in range(M)] for X in range(M) for Y in range(M)])
        Rcol = flint.fmpz_mat([[T[X][Y][U][V] for Y in range(M) for V in range(M)] for X in range(M) for U in range(M)])
        rr, rc = Rrow.rank(), Rcol.rank()
        brow = sum(p(s + i) for i in range(2 * h + 1) if s + i >= 0)
        bcol = sum(p(2 * h - s - i) for i in range(2 * h + 1) if 2 * h - s - i >= 0)
        rows_out.append((s, rr, brow, rc, bcol))
        total_min += min(rr, rc)
    for t in rows_out: print("  s=%3d rank_row=%5d (bound %5d)  rank_col=%5d (bound %5d)" % t)
    print(f"  => P({M}) <= sum_s min(rank_row,rank_col) = {total_min};  N^2 = {M*M};  "
          f"log_N = {__import__('math').log(total_min)/__import__('math').log(M):.3f}")
    return total_min

if __name__ == "__main__":
    # cross-check recursion (O2) against definition (47)
    import itertools
    bad = 0
    for parts in [(1,), (2,), (1, 1), (2, 1), (3,), (2, 2), (3, 1), (1, 1, 1)]:
        for E in range(0, 5):
            for sigma in parts_list(E):
                for r in range(-6, 7):
                    d1 = dict(apply_mono(parts, r, sigma)); d2 = apply_direct(parts, r, sigma)
                    if d1 != d2: bad += 1
    print("O2 recursion vs (47) mismatches:", bad)
    for h in map(int, sys.argv[1:]):
        run(h)
