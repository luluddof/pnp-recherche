"""Independent re-implementation (agent G5b) of Appendix B pieces, written from the paper's
definitions (47),(48),(53), Lemma 10; does not use heis.py.
Operator t^mu[r]: sum over nonzero k_i (one per part), sum k_i = -r, coefficient prod C(-k_i-1, mu_i-1),
differentiations nabla_k = k d/dt_k (k>0) first, then multiplications t_{-k} (k<0).
"""
import sys, itertools
from functools import lru_cache
from math import comb, factorial
from collections import defaultdict
import flint

def gbin(n, k):  # generalized binomial C(n,k), n any integer, k>=0
    if k < 0: return 0
    num = 1
    for t in range(k): num *= (n - t)
    return num // factorial(k)

@lru_cache(None)
def parts(n):
    if n < 0: return ()
    out = []
    def rec(rem, mx, cur):
        if rem == 0: out.append(tuple(cur)); return
        for k in range(min(rem, mx), 0, -1): rec(rem - k, k, cur + [k])
    rec(n, n, [])
    return tuple(out)

@lru_cache(None)
def pidx(n): return {m: i for i, m in enumerate(parts(n))}

def srt(t): return tuple(sorted(t, reverse=True))

@lru_cache(None)
def multpoly(negparts, L):
    """sum over l_i >= 1 with sum l_i = L of prod C(l_i - 1, mu_i - 1) * t^{(l_i)}"""
    out = defaultdict(int)
    m = len(negparts)
    if m == 0:
        return ((( ), 1),) if L == 0 else ()
    def rec(j, rem, coef, mono):
        if j == m - 1:
            l = rem
            if l < 1: return
            c = coef * gbin(l - 1, negparts[j] - 1)
            if c: out[srt(mono + (l,))] += c
            return
        for l in range(1, rem - (m - 1 - j) + 1):
            c = coef * gbin(l - 1, negparts[j] - 1)
            if c: rec(j + 1, rem - l, c, mono + (l,))
    rec(0, L, 1, ())
    return tuple((k, v) for k, v in out.items() if v)

@lru_cache(None)
def apply(mu, r, sigma):
    """t^mu [r] applied to t^sigma; returns dict mono->coef (as tuple of items)"""
    out = defaultdict(int)
    def rec(j, f, neg, sumpos):
        if j == len(mu):
            L = r + sumpos
            G = multpoly(tuple(sorted(neg)), L)
            for mono, c in f.items():
                for g, d in G:
                    out[srt(mono + g)] += c * d
            return
        n = mu[j]
        # negative choice
        rec(j + 1, f, neg + [n], sumpos)
        # positive choices k
        maxw = max((sum(m) for m in f), default=0)
        for k in range(1, maxw + 1):
            c0 = gbin(-k - 1, n - 1)
            g = defaultdict(int)
            for mono, c in f.items():
                cnt = mono.count(k)
                if cnt:
                    l = list(mono); l.remove(k)
                    g[tuple(l)] += c * c0 * k * cnt
            g = {a: b for a, b in g.items() if b}
            if g: rec(j + 1, g, neg, sumpos + k)
    rec(0, {sigma: 1}, [], 0)
    return tuple((k, v) for k, v in out.items() if v)

def apply_poly(poly, r, vec):
    out = defaultdict(int)
    for mono, c in poly.items():
        for sig, d in vec.items():
            for k, v in apply(mono, r, sig):
                out[k] += c * d * v
    return {k: v for k, v in out.items() if v}

def pmul(f, g):
    out = defaultdict(int)
    for a, x in f.items():
        for b, y in g.items(): out[srt(a + b)] += x * y
    return dict(out)

def zee(mu):
    z = 1
    for k in set(mu): z *= k ** mu.count(k) * factorial(mu.count(k))
    return z

def run(h):
    P = parts(h); M = len(P)
    a = {}
    for mu in P:
        f = {(): 1}
        for n in mu:
            f = pmul(f, {(j,): (-1) ** (n - 1) * comb(n, j) for j in range(1, n + 1)})
        a[mu] = {k: v for k, v in f.items() if v}
    b = {nu: {nu: factorial(h) // zee(nu)} for nu in P}
    # Lemma 10 check
    l10 = all(sum(c * dict(apply_poly({m: 1}, -h, {nu: factorial(h)//zee(nu)})).get((), 0)
                  for m, c in a[mu].items()) == factorial(h) * (mu == nu) for mu in P for nu in P)
    S = list(range(-2 * h, 2 * h + 1))
    A = {s: [[[[0]*M for _ in range(M)] for _ in range(M)] for _ in range(M)] for s in S}  # A[s][X][Y][U][V]
    for s in S:
        for i in range(0, 2 * h + 1):
            r = h - s - i
            if h + r < 0: continue
            sg = (-1) ** i
            for vi, V in enumerate(P):
                mids = [apply_poly(b[V], r, {Y: 1}) for Y in P]
                for xi, X in enumerate(P):
                    ai = {m: c * comb(sum(m) + h, i) for m, c in a[X].items()}
                    ai = {m: c for m, c in ai.items() if c}
                    for yi in range(M):
                        res = apply_poly(ai, -r, mids[yi])
                        for mono, c in res.items():
                            A[s][xi][yi][pidx(h)[mono]][vi] += sg * c
    hf = factorial(h)
    ok = all(sum(A[s][x][y][u][v] for s in S) == hf * (x == v) * (y == u)
             for x in range(M) for y in range(M) for u in range(M) for v in range(M))
    # also check that with b <-> a order exchanged (Lemma 9a) at random spot? (compute col factor directly)
    tot = 0; lines = []
    for s in S:
        T = A[s]
        rr = flint.fmpz_mat([[T[x][y][u][v] for u in range(M) for v in range(M)] for x in range(M) for y in range(M)]).rank()
        rc = flint.fmpz_mat([[T[x][y][u][v] for y in range(M) for v in range(M)] for x in range(M) for u in range(M)]).rank()
        brow = sum(len(parts(s + i)) for i in range(2 * h + 1))
        bcol = sum(len(parts(2 * h - s - i)) for i in range(2 * h + 1))
        lines.append((s, rr, brow, rc, bcol, rr <= brow and rc <= bcol))
        tot += min(rr, rc)
    print(f"h={h} M={M} Lemma10={l10} identity={ok} sum_min={tot} all_ranks_within_bounds={all(l[-1] for l in lines)}")
    for l in lines: print("   s=%3d row %4d<=%5d col %4d<=%5d %s" % l)
    nz = [s for s in S if any(A[s][x][y][u][v] for x in range(M) for y in range(M) for u in range(M) for v in range(M))]
    print("   nonzero pieces s:", nz)
    sys.stdout.flush()
    return A, P

if __name__ == "__main__":
    for h in map(int, sys.argv[1:]): run(h)
