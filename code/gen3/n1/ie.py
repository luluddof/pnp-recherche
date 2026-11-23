"""N1 brick: inclusion-exclusion (IE) over marked clauses for #3-SAT.
RECONSTRUCTION from a web snippet only; not the paper's algorithm.

Identity:  #SAT(F) = sum_{S subset clauses} (-1)^|S| N(S),
N(S) = #assignments violating every clause of S = 2^{n-|vars(S)|} if consistent else 0.
Equivalent "moment" form: let r(x) = #clauses violated by x, a_r = #{x : r(x)=r}.
  b_k := sum_{|S|=k} N(S) = sum_r a_r C(r,k)
  T_t := sum_{k<=t} (-1)^k b_k = a_0 + (-1)^t * sum_{r>t} a_r C(r-1,t)     (*)
so truncation at t is EXACT on every assignment with r(x) <= t, and the whole
error comes from assignments violating > t clauses.
"""
import random, itertools, math
from math import comb
import numpy as np


def rand3sat(n, m, rng):
    cls = []
    for _ in range(m):
        vs = rng.sample(range(1, n + 1), 3)
        cls.append([v if rng.random() < .5 else -v for v in vs])
    return cls


def viol_counts(n, clauses):
    """r(x) for all 2^n assignments (bit v-1 of x = value of var v)."""
    N = 1 << n
    x = np.arange(N, dtype=np.int64)
    r = np.zeros(N, dtype=np.int16)
    bits = [((x >> (v - 1)) & 1).astype(bool) for v in range(1, n + 1)]
    for c in clauses:
        viol = np.ones(N, dtype=bool)
        for l in c:
            b = bits[abs(l) - 1]
            viol &= (~b) if l > 0 else b   # literal false
        r += viol
    return r


def hist(n, clauses):
    r = viol_counts(n, clauses)
    return np.bincount(r, minlength=len(clauses) + 1).astype(object)


def T_from_hist(a, t):
    return sum((-1) ** k * sum(int(a[r]) * comb(r, k) for r in range(len(a))) for k in range(t + 1))


def T_direct(n, clauses, t):
    """Truncated IE by explicit enumeration of consistent clause sets S, |S|<=t."""
    m = len(clauses)
    total = 0
    def rec(i, k, fixed):
        nonlocal total
        total += (-1) ** k * 2 ** (n - len(fixed))
        if k == t:
            return
        for j in range(i, m):
            nf = dict(fixed); ok = True
            for l in clauses[j]:          # violating clause j forces each literal false
                v, val = abs(l), (l < 0)
                if nf.get(v, val) != val: ok = False; break
                nf[v] = val
            if ok:
                rec(j + 1, k + 1, nf)
    rec(0, 0, {})
    return total


def count_consistent(n, clauses, t):
    m = len(clauses); cnt = [0] * (t + 1)
    def rec(i, k, fixed):
        cnt[k] += 1
        if k == t: return
        for j in range(i, m):
            nf = dict(fixed); ok = True
            for l in clauses[j]:
                v, val = abs(l), (l < 0)
                if nf.get(v, val) != val: ok = False; break
                nf[v] = val
            if ok: rec(j + 1, k + 1, nf)
    rec(0, 0, {})
    return cnt


def rmax_maxsat(n, clauses):
    """max_x r(x) exactly via MaxSAT (RC2): y_j -> clause j violated."""
    from pysat.examples.rc2 import RC2
    from pysat.formula import WCNF
    w = WCNF()
    for j, c in enumerate(clauses):
        y = n + 1 + j
        for l in c:
            w.append([-y, -l])
        w.append([y], weight=1)
    with RC2(w) as s:
        s.compute()
        return len(clauses) - s.cost


def is_sat(n, clauses):
    from pysat.solvers import Minisat22
    with Minisat22(bootstrap_with=clauses) as s:
        return s.solve()


def t_decide(a):
    """Smallest t at which Bonferroni certifies the answer:
       even t: T_t >= a0, so T_t < 1 => UNSAT;  odd t: T_t <= a0, so T_t >= 1 => SAT."""
    m = len(a) - 1
    for t in range(m + 1):
        T = T_from_hist(a, t)
        if t % 2 == 0 and T < 1: return t, 'UNSAT'
        if t % 2 == 1 and T >= 1: return t, 'SAT'
    return m, ('SAT' if a[0] else 'UNSAT')
