"""G9a: Appendix-B pieces with index set restricted to partitions with parts in S.
Mode 'full'  : operators of the full Heisenberg VOA, x,y,u,v restricted to Par_S(h).
               (valid decomposition by restriction; we measure the EXACT ranks.)
Mode 'comp'  : intermediate states projected onto W_S = Q[t_k : k in S] (compressed
               operators pi_S a[r] pi_S) -- tests whether Lemma 9 survives on the sub-Fock space.
"""
import sys, math
from functools import lru_cache
from math import comb, factorial
from collections import defaultdict
import flint

def partitions(n, maxp=None):
    if maxp is None: maxp = n
    if n == 0:
        yield (); return
    for k in range(min(n, maxp), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest

def gbinom(k, j):
    if j < 0: return 0
    num = 1
    for t in range(j): num *= (k - t)
    return num // factorial(j)

def add_part(sigma, l): return tuple(sorted(sigma + (l,), reverse=True))
def rem_part(sigma, k):
    s = list(sigma); s.remove(k); return tuple(s)

@lru_cache(None)
def apply_mono(parts, r, sigma):
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

def apply_poly(poly, r, vec, S=None):
    """poly[r] applied to vector vec (dict mono->coef); optional projection onto W_S."""
    out = defaultdict(int)
    for sigma, c0 in vec.items():
        for mono, c in poly.items():
            for m, v in apply_mono(mono, r, sigma):
                out[m] += c0 * c * v
    if S is not None:
        return {m: v for m, v in out.items() if v and all(q in S for q in m)}
    return {m: v for m, v in out.items() if v}

def poly_mul(f, g):
    out = defaultdict(int)
    for a, x in f.items():
        for b, y in g.items():
            out[tuple(sorted(a + b, reverse=True))] += x * y
    return {k: v for k, v in out.items() if v}

def zee(mu):
    z = 1
    for k in set(mu):
        mk = mu.count(k); z *= k ** mk * factorial(mk)
    return z

def run(h, S, mode):
    P = [mu for mu in partitions(h) if S is None or all(q in S for q in mu)]
    M = len(P)
    y = {n: {(j,): comb(n, j) for j in range(1, n + 1)} for n in range(1, h + 1)}
    a = {}
    for mu in P:
        f = {(): 1}
        for part in mu:
            f = poly_mul(f, {k: (-1) ** (part - 1) * v for k, v in y[part].items()})
        a[mu] = f
    b = {nu: {nu: factorial(h) // zee(nu)} for nu in P}
    def a_i(mu, i):
        return {k: v * comb(sum(k) + h, i) for k, v in a[mu].items() if comb(sum(k) + h, i)}
    proj = S if mode == 'comp' else None
    Srange = range(-2 * h, 2 * h + 1)
    idx = {mu: i for i, mu in enumerate(P)}
    A = {}; agree = True
    for s in Srange:
        T = {}
        for X in P:
            ais = [a_i(X, i) for i in range(2 * h + 1)]
            for V in P:
                for Y in P:
                    vy = {Y: 1}
                    orig = defaultdict(int); exch = defaultdict(int)
                    for i in range(2 * h + 1):
                        sg = (-1) ** i
                        w1 = 2 * h - s - i; w2 = s + i
                        do_orig = w1 >= 0 and (mode == 'comp' or s >= 0)
                        do_exch = w2 >= 0 and (mode == 'comp' or s < 0)
                        if do_orig:
                            t = apply_poly(b[V], h - s - i, vy, proj)
                            t = apply_poly(ais[i], s - h + i, t, None)
                            for m, v in t.items():
                                if m in idx: orig[m] += sg * v
                        if do_exch:
                            t = apply_poly(ais[i], s - h + i, vy, proj)
                            t = apply_poly(b[V], h - s - i, t, None)
                            for m, v in t.items():
                                if m in idx: (exch if mode == 'comp' else orig)[m] += sg * v
                    if mode == 'comp':
                        for U in P:
                            if orig.get(U, 0) != exch.get(U, 0): agree = False
                    for U in P:
                        T[(idx[X], idx[Y], idx[U], idx[V])] = (orig.get(U, 0), exch.get(U, 0))
        A[s] = T
    hf = factorial(h)
    def ident(which):
        ok = True
        for X in range(M):
            for Y in range(M):
                for U in range(M):
                    for V in range(M):
                        tot = sum(A[s][(X, Y, U, V)][which] for s in Srange)
                        if tot != hf * (X == V) * (Y == U): ok = False
        return ok
    name = 'all' if S is None else sorted(q for q in S if q <= h)
    if mode == 'comp':
        print(f"[comp] h={h} S={name} N={M}: orders agree (Lemma 9a): {agree};  sum_s A_s=h!*crossing: b-first {ident(0)}, a-first {ident(1)}")
        return
    ok = ident(0)
    tot_min = 0; sum_rr = 0; maxcost = 0.0; maxcost_pos = 0.0
    for s in Srange:
        T = A[s]
        Rrow = flint.fmpz_mat([[T[(X, Y, U, V)][0] for U in range(M) for V in range(M)] for X in range(M) for Y in range(M)])
        Rcol = flint.fmpz_mat([[T[(X, Y, U, V)][0] for Y in range(M) for V in range(M)] for X in range(M) for U in range(M)])
        rr, rc = Rrow.rank(), Rcol.rank()
        tot_min += min(rr, rc); sum_rr += rr
        if rr and rc: maxcost = max(maxcost, math.log2(rr) + math.log2(rc))
    lam = math.log2(M) if M > 1 else float('nan')
    print(f"[full] h={h} S={name} N={M} ident={ok} sum_min_rank={tot_min} (N^2={M*M}) sum_rank_row={sum_rr} "
          f"log_N(sum_min)={math.log(tot_min)/math.log(M) if M>1 else float('nan'):.3f} rho_eff=max_s(c_row+c_col)/(4 lam)={maxcost/(4*lam):.3f}")

if __name__ == "__main__":
    mode = sys.argv[1]; Sname = sys.argv[2]
    hs = list(map(int, sys.argv[3:]))
    for h in hs:
        S = {'all': None, 'pow2': {1, 2, 4, 8, 16}, 'sq': {1, 4, 9, 16}, 'le2': {1, 2}, 'le3': {1, 2, 3},
             'odd': set(range(1, 40, 2)), 'ge2': set(range(2, 40)), 'even': set(range(2, 40, 2)), '13': {1, 3}}[Sname]
        run(h, S, mode)
