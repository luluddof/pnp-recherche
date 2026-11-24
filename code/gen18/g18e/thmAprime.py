#!/usr/bin/env python3
"""G18e: Theorem A' (unbalanced-group flattening).  Exact arithmetic only.
 (1) sanity: rank of restricted direction terms on a toy instance, exact over Q (sympy), vs 2^{c_i};
 (2) product labels {0,1}^u (halves of [r] merged), uniform ordered-pair loads: bound for threshold cuts,
     compared with Theorem A (2M/4^u);  (3) brute force over ALL cuts for u = 2, 3."""
import itertools, math, random
from fractions import Fraction as Fr
import sympy as sp

def bound_labels(u, A, B):
    """labels in {0,1}^u as ints; one wire per ORDERED pair of labels (M = 4^u). returns (w, max_i c_i)."""
    w = 2 * len(A) * len(B); cmax = 0
    for i in range(u):
        c = 0
        for g in (0, 1):
            e1 = 2 * len(B) * sum(1 for x in A if (x >> i) & 1 == g)
            e2 = 2 * len(A) * sum(1 for x in B if (x >> i) & 1 == g)
            c += min(e1, e2)
        cmax = max(cmax, c)
    return w, cmax

print("(2) threshold cuts A={s>=t1}, B={s<=t0}, s = Hamming weight; log2 N >= frac * M")
print("  u | ThmA 2/4^u | ThmA' best threshold (t0,t1) | majority formula C(u-1,(u-1)/2)/2^u | A'*sqrt(u)")
for u in list(range(1, 15)) + [21, 31, 41, 61, 101, 201, 401]:
    best = (Fr(0), None)
    cnt = [math.comb(u, s) for s in range(u + 1)]
    for t0 in range(0, u):
        for t1 in range(t0 + 1, u + 1):
            a = sum(cnt[s] for s in range(t1, u + 1)); b = sum(cnt[s] for s in range(0, t0 + 1))
            # number of labels in A with bit i = 1: sum_s C(u-1,s-1)
            a1 = sum(math.comb(u - 1, s - 1) for s in range(max(t1, 1), u + 1)); a0 = a - a1
            b1 = sum(math.comb(u - 1, s - 1) for s in range(1, t0 + 1)); b0 = b - b1
            w = 2 * a * b; c = min(2 * b * a1, 2 * a * b1) + min(2 * b * a0, 2 * a * b0)
            val = Fr(w - c, 4 ** u)
            if val > best[0]: best = (val, (t0, t1))
    maj = Fr(math.comb(u - 1, (u - 1) // 2), 2 ** u) if u % 2 else None
    print(f"  {u:3d} | {2/4**u:.3e} | {float(best[0]):.5f} {best[1]} | {float(maj) if maj else float('nan'):.5f} | {float(best[0])*math.sqrt(u):.4f}")
print("  limit of A'*sqrt(u) for the majority cut: 1/sqrt(2 pi) =", 1 / math.sqrt(2 * math.pi))

print("(3) brute force over all disjoint (A,B) of labels, uniform ordered pairs")
for u in (1, 2, 3):
    n = 2 ** u; best = (-1, None)
    for assign in itertools.product((0, 1, 2), repeat=n):
        A = [x for x in range(n) if assign[x] == 1]; B = [x for x in range(n) if assign[x] == 2]
        if not A or not B: continue
        w, c = bound_labels(u, A, B)
        if w - c > best[0]: best = (w - c, (A, B, w, c))
    print(f"  u={u}: max (w - c) = {best[0]} of M = {4**u}  ratio {best[0]/4**u:.5f}   cut {best[1]}   [Thm A: {2} of {4**u}]")

print("(1) sanity on a toy: labels {0,1}^2, one wire per unordered pair of distinct labels (M=6)")
labels = [0, 1, 2, 3]; wires = list(itertools.combinations(labels, 2))
random.seed(1)
def restricted_term_rank(i, C1, C2):
    cut = [(a, b) if a in C1 else (b, a) for (a, b) in wires if (a in C1 and b in C2) or (a in C2 and b in C1)]
    w = len(cut)
    # term = prod_g tau_g; tau_g random integer function of the cut-wire end bits located in g
    rows = list(itertools.product((0, 1), repeat=w)); Mx = sp.zeros(2 ** w, 2 ** w)
    taus = {}
    for g in (0, 1):
        r_idx = [k for k, (a, b) in enumerate(cut) if (a >> i) & 1 == g]
        c_idx = [k for k, (a, b) in enumerate(cut) if (b >> i) & 1 == g]
        taus[g] = (r_idx, c_idx, {})
    for ri, R in enumerate(rows):
        for ci, Cc in enumerate(rows):
            v = 1
            for g in (0, 1):
                r_idx, c_idx, tab = taus[g]
                key = (tuple(R[k] for k in r_idx), tuple(Cc[k] for k in c_idx))
                if key not in tab: tab[key] = random.randint(-5, 5)
                v *= tab[key]
            Mx[ri, ci] = v
    c = sum(min(len(taus[g][0]), len(taus[g][1])) for g in (0, 1))
    return w, c, Mx.rank()
ok = True; bestb = 0
for assign in itertools.product((0, 1, 2), repeat=4):
    C1 = {x for x in labels if assign[x] == 1}; C2 = {x for x in labels if assign[x] == 2}
    if not C1 or not C2: continue
    cs = []
    for i in (0, 1):
        w, c, rk = restricted_term_rank(i, C1, C2); ok &= (rk <= 2 ** c); cs.append((c, rk))
    bestb = max(bestb, w - max(c for c, _ in cs))
print("   all random restricted terms have rank <= 2^{c_i}:", ok, "| best A' bound on this toy: log2 N >=", bestb, "(Thm A: 1; trivial upper bound: 4)")
