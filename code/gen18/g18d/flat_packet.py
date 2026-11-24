# (1) Smith profile robustness (other t, zeta).  (2) Flattening of ONE packet tensor across the two groups of
# coordinate i, 2-adic Smith invariants, vs number of mode monomials of weight < p.
import sys, itertools
from fractions import Fraction
sys.path.insert(0, '.')
import importlib.util
K = 1200; MOD = 1 << K
def v2(n):
    if n == 0: return K
    return (n & -n).bit_length() - 1
def to_mod(fr):
    assert fr.denominator % 2 == 1
    return (fr.numerator * pow(fr.denominator, -1, MOD)) % MOD
def smith_vals(A, cap):
    A = [row[:] for row in A]; rows = list(range(len(A))); cols = list(range(len(A[0]))); vals = []
    while rows and cols:
        best = None
        for r in rows:
            Ar = A[r]
            for c in cols:
                e = Ar[c]
                if e:
                    v = (e & -e).bit_length() - 1
                    if best is None or v < best[0]:
                        best = (v, r, c)
                        if v == (vals[-1] if vals else 0): break
            if best and best[0] == (vals[-1] if vals else 0): break
        if best is None or best[0] >= cap: break
        v, pr, pc = best; vals.append(v)
        inv = pow(A[pr][pc] >> v, -1, MOD); P = A[pr]
        for r in rows:
            if r == pr: continue
            e = A[r][pc]
            if e == 0: continue
            m = ((e >> v) * inv) % MOD; Ar = A[r]
            for c in cols:
                if P[c]: Ar[c] = (Ar[c] - m * P[c]) % MOD
        rows.remove(pr); cols.remove(pc)
    return vals
def Wmat(a, T):
    L = 1 << a; S = range(-L, L)
    Cp0 = 1
    for q in S:
        if q: Cp0 *= -q
    CT = 1
    for q in S: CT *= (T - q)
    return [[to_mod(Fraction(CT, Cp0 * (T - (y - x)))) for y in range(L)] for x in range(L)]
print("## (1) crossing-block Smith valuations, other nodes / zeta")
for (a, t, zeta) in [(5, 3, 1), (5, -1, 1), (5, 6, 3), (5, 2, 5), (5, 12, 7), (4, 1, 2), (4, 1, 4), (5, -4, 1)]:
    L = 1 << a; i = v2(t); T = t + (1 << (a + 1)) * zeta; W = Wmat(a, T)
    xs = [x for x in range(L) if (x >> i) & 1 == 0]
    vals = sorted(smith_vals([[W[x][y] for y in xs] for x in xs], 400))
    pred = sorted([(a - i + 1) + v2(zeta) + 2 * k + 2 * sum(v2(m) for m in range(1, k + 1)) for k in range(L >> (i + 1))] * (1 << i))
    print(f"a={a} t={t} zeta={zeta} match_formula(+v2(zeta))={vals == pred} vals[:8]={vals[:8]}")
print("## (2) flattening of one packet tensor across groups; #invariants < p  vs  #mode monomials of weight < p")
def count_monomials(weights, p):
    # number of exponent tuples with sum w*beta < p (one variable per entry of weights)
    dp = [0] * p; dp[0] = 1
    for w in weights:
        for s in range(w, p): dp[s] += dp[s - w]
    return sum(dp)
for (a, i) in [(2, 0), (2, 1), (3, 0), (3, 1), (3, 2)]:
    L = 1 << a; t = 1 << i; T = t + (1 << (a + 1)) * 1; W = Wmat(a, T)
    rg = [x for x in range(L) if (x >> i) & 1 == 0]; cg = [y for y in range(L) if (y >> i) & 1 == 1]   # group g
    rh = [x for x in range(L) if (x >> i) & 1 == 1]; ch = [y for y in range(L) if (y >> i) & 1 == 0]   # group g'
    h = L // 2
    def bits(m, n): return [(m >> j) & 1 for j in range(n)]
    rowidx = list(itertools.product(range(1 << h), repeat=2))   # (Y_g, X_g)
    colidx = rowidx
    A = []
    for (Yg, Xg) in rowidx:
        yg = bits(Yg, h); xg = bits(Xg, h); row = []
        for (Yh, Xh) in colidx:
            yh = bits(Yh, h); xh = bits(Xh, h)
            X = [0] * L; Y = [0] * L
            for k, c in enumerate(cg): X[c] = xg[k]
            for k, c in enumerate(ch): X[c] = xh[k]
            for k, r in enumerate(rg): Y[r] = yg[k]
            for k, r in enumerate(rh): Y[r] = yh[k]
            val = 1
            for x in range(L):
                s = 0
                for y in range(L):
                    if X[y]: s += W[x][y]
                val = (val * (s if Y[x] else 1 - s)) % MOD
                if val == 0: break
            row.append(val)
        A.append(row)
    cap = 40
    vals = sorted(smith_vals(A, cap))
    D = L >> i
    sm = [(a - i + 1) + 2 * k + 2 * sum(v2(m) for m in range(1, k + 1)) for k in range(D // 2)] * (1 << i) * 2  # Smith modes, both directions
    pm = [(a - i + 1) + b for b in range(cap)] * (1 << i) * 2      # paper's modes (true weights), both directions
    out = []
    for p in (2, 4, 6, 8, 10, 12, 16, 20, 24, 32, 40):
        out.append(f"p={p}:{sum(1 for v in vals if v < p)}/{count_monomials(sm, p)}/{count_monomials(pm, p)}")
    print(f"a={a} i={i} L={L} size={len(A)}  (flattening rank mod 2^p)/(Smith-mode monomials)/(paper-mode monomials): {' '.join(out)}")
    print(f"     first invariants: {vals[:40]}")
    sys.stdout.flush()
