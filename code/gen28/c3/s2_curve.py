# Theorem S'' bound vs singlet family (n=2) at finite k.  eta = log2(q)/k (N=2^k), delta = -log2(p/N^2)/k.
from math import comb, log2
from fractions import Fraction
def err_bound(k, d):
    # exact rational: p/N^2 >= sum_m P(m) P(Bin(m,1/2) > (m+d)/2), m ~ Bin(k,1/2)
    s = Fraction(0)
    for m in range(k+1):
        tail = sum(comb(m, i) for i in range((m+d)//2 + 1, m+1))
        s += Fraction(comb(k, m), 2**k) * Fraction(tail, 2**m)
    return s
def singlet_best(k, q):
    # windows [a, a+d] of the singlet number, cost q(d) <= sum_{u<=d} C(k,u)3^u (upper bound, exact for these windows generically)
    best = None
    pm = [Fraction(comb(k, m) * 3**(k-m), 4**k) for m in range(k+1)]
    for d in range(k+1):
        qd = sum(comb(k, u) * 3**u for u in range(d+1))
        if qd > q: break
        e = min(1 - sum(pm[a:a+d+1]) for a in range(k+1))
        best = e
    return best
import sys
for k in [20, 40, 80, 160]:
    print(f"k={k}")
    for d in range(1, k//3, max(1, k//20)):
        q = 3**d  # weakest cost compatible with degree d
        eb = err_bound(k, d)
        eta = d*log2(3)/k
        db = -log2(float(eb))/k if eb > 0 else float('inf')
        es = singlet_best(k, 2**(eta*k)) if k <= 80 else None
        ds = (-log2(float(es))/k if es else 0.0) if es is not None else float('nan')
        print(f"  d={d:3d} eta={eta:.3f}  delta_bound={db:.4f} (/eta^2={db/eta**2:.3f})   singlet delta at same q={ds:.4f}")
