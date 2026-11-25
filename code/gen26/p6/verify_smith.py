# Independent check of the Smith-profile theorem for the crossing block of the Cauchy kernel (D.7).
# E[x][y] = C(T)/(C'(0)(T-(y-x))), x,y in [0,L) with bit_i(x)=bit_i(y)=b, T = t + 2^{a+1} zeta, v2(t)=i.
# Claim: Smith valuations over Z_(2) = {(a-i+1)+v2(zeta)+2k+2v2(k!) : 0<=k<2^{a-i-1}}, each with multiplicity 2^i.
import sys, itertools
from fractions import Fraction
K = 900; MOD = 1 << K
def v2(n):
    n = abs(n)
    return K if n == 0 else (n & -n).bit_length() - 1
def smith_vals(A):
    A = [r[:] for r in A]; rows = list(range(len(A))); cols = list(range(len(A[0]))); out = []
    while rows and cols:
        best = None
        for r in rows:
            for c in cols:
                e = A[r][c]
                if e:
                    v = (e & -e).bit_length() - 1
                    if best is None or v < best[0]: best = (v, r, c)
        if best is None or best[0] >= K - 100: break
        v, pr, pc = best; out.append(v); inv = pow(A[pr][pc] >> v, -1, MOD)
        for r in rows:
            if r != pr and A[r][pc]:
                m = ((A[r][pc] >> v) * inv) % MOD
                A[r] = [(A[r][c] - m * A[pr][c]) % MOD for c in range(len(A[r]))]
        rows.remove(pr); cols.remove(pc)
    return sorted(out)
def block(a, t, zeta, b):
    L = 1 << a; S = range(-L, L); T = t + (1 << (a + 1)) * zeta
    Cp0 = 1
    for q in S:
        if q: Cp0 *= -q
    CT = 1
    for q in S: CT *= (T - q)
    i = v2(t); B = [x for x in range(L) if (x >> i) & 1 == b]
    def tm(fr): return (fr.numerator * pow(fr.denominator, -1, MOD)) % MOD
    return [[tm(Fraction(CT, Cp0 * (T - (y - x)))) for y in B] for x in B], CT, Cp0
def formula(a, i, zeta):
    return sorted([(a - i + 1) + v2(zeta) + 2 * k + 2 * sum(v2(m) for m in range(1, k + 1))
                   for k in range(1 << (a - i - 1))] * (1 << i))
ok = 0; bad = 0
for a in range(1, 8):
    L = 1 << a
    ts = [t for t in range(-L, L) if t != 0 and v2(t) < a]
    if a >= 6: ts = [t for t in ts if abs(t) <= 9 or t in (L - 1, L // 2, -L // 2, 3 * L // 8)]
    for t in ts:
        for zeta in ((1, 2, 3, -5, 12) if a <= 4 else (1, 6)):
            for b in ((0, 1) if a <= 5 else (0,)):
                E, CT, Cp0 = block(a, t, zeta, b)
                assert v2(CT) - v2(Cp0) == a + 1 + v2(zeta)
                if smith_vals(E) == formula(a, v2(t), zeta): ok += 1
                else:
                    bad += 1; print("MISMATCH", a, t, zeta, b, smith_vals(E)[:10], formula(a, v2(t), zeta)[:10])
    print(f"a={a} done: ok={ok} bad={bad}"); sys.stdout.flush()
# brute force: minimal minor valuation achieved by prefix sets, small case
def minor_val(a, t, zeta, R, C):
    T = t + (1 << (a + 1)) * zeta; i = v2(t)
    return sum(v2(x - y) for x, y in itertools.combinations(R, 2)) + sum(v2(x - y) for x, y in itertools.combinations(C, 2)) \
        - sum(v2(T + x - y) for x in R for y in C)
for (a, t) in [(3, 1), (3, 2), (3, 4), (4, 2), (4, 4), (4, 6)]:
    i = v2(t); B = [x for x in range(1 << a) if (x >> i) & 1 == 0]
    for k in range(1, len(B) + 1):
        if len(B) > 8 and k > 4: break
        m = min(minor_val(a, t, 1, R, C) for R in itertools.combinations(B, k) for C in itertools.combinations(B, k))
        f = formula(a, i, 1); pred = sum(f[:k]) - k * (a + 1)
        assert m == pred, (a, t, k, m, pred)
    print("brute-force minors ok", a, t)
