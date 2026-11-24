# 2-adic Smith profile of the CROSSING block of App. D's Cauchy kernel W(t + alpha*zeta), one full packet.
# Rows x, columns y with bit_i(x) == bit_i(y) (different groups in coordinate i = v2(t)).
# Compares #invariants of valuation <= w with the number of modes of true weight <= w (b+1+a-i per tile).
import sys
from fractions import Fraction
K = 3000
MOD = 1 << K
def v2(n):
    if n == 0: return K
    return (n & -n).bit_length() - 1
def to_mod(fr):
    num, den = fr.numerator, fr.denominator
    assert den % 2 == 1
    return (num * pow(den, -1, MOD)) % MOD
def smith_vals(A):
    A = [row[:] for row in A]; R = len(A); Cn = len(A[0]); vals = []
    rows = list(range(R)); cols = list(range(Cn))
    while rows and cols:
        best = None
        for r in rows:
            for c in cols:
                v = v2(A[r][c])
                if best is None or v < best[0]: best = (v, r, c)
        v, pr, pc = best
        if v >= K - 200: break
        vals.append(v)
        piv = A[pr][pc] >> v
        inv = pow(piv, -1, MOD)
        for r in rows:
            if r == pr: continue
            e = A[r][pc]
            if e == 0: continue
            m = ((e >> v) * inv) % MOD
            for c in cols:
                A[r][c] = (A[r][c] - m * A[pr][c]) % MOD
        rows.remove(pr); cols.remove(pc)
    return vals
def block(a, t, zeta, bit):
    L = 1 << a; i = v2(t); alpha = 1 << (a + 1)
    S = list(range(-L, L)); T = t + alpha * zeta
    Cp0 = 1
    for q in S:
        if q != 0: Cp0 *= (0 - q)
    CT = 1
    for q in S: CT *= (T - q)
    xs = [x for x in range(L) if (x >> i) & 1 == bit]
    ys = [y for y in range(L) if (y >> i) & 1 == bit]
    return [[to_mod(Fraction(CT, Cp0 * (T - (y - x)))) for y in ys] for x in xs]
for a in (4, 5, 6):
    L = 1 << a
    for i in range(0, a):
        t = 1 << i
        for zeta in (1,):
            vals = sorted(smith_vals(block(a, t, zeta, 0)))
            D = L >> i; ntile = 1 << i     # tiles in this block
            # modes of true weight b+1+a-i, b = 0..D/2-1, per tile
            prof = []
            for w in (a - i + 1, a - i + 2, a - i + 4, a - i + 8, a - i + 16):
                Rt = sum(1 for v in vals if v <= w)
                Rm = min(ntile * max(0, w - (a - i)), L // 2)
                prof.append(f"w={w}:{Rt}/{Rm}")
            # per-tile slope: valuation of k-th invariant ~ slope*k/ntile
            print(f"a={a} i={i} D={D} tiles={ntile} n_inv={len(vals)} min={vals[0]} max={vals[-1]} "
                  f"(true rank<=w)/(modes<=w): {' '.join(prof)}  vals[:12]={vals[:12]}")
            sys.stdout.flush()
