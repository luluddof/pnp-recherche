# C1 test: flattening F of one packet tensor across groups (g | g') vs tensor product of the two one-direction
# "Boolean product" matrices M[S, X'] = prod_{x in S} (E X')_x  (S subset of target rows, X' in {0,1}^{source cols}).
import sys, itertools
from fractions import Fraction
K = 700; MOD = 1 << K
def v2(n):
    n = abs(n); return K if n == 0 else (n & -n).bit_length() - 1
def smith_vals(A, cap):
    A = [r[:] for r in A]; rows = set(range(len(A))); cols = list(range(len(A[0]))); out = []
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
            if best and out and best[0] == out[-1]: break
        if best is None or best[0] >= cap: break
        v, pr, pc = best; out.append(v); inv = pow(A[pr][pc] >> v, -1, MOD); P = A[pr]
        nz = [c for c in cols if P[c]]
        for r in rows:
            if r != pr and A[r][pc]:
                m = ((A[r][pc] >> v) * inv) % MOD; Ar = A[r]
                for c in nz: Ar[c] = (Ar[c] - m * P[c]) % MOD
        rows.discard(pr); cols.remove(pc)
    return sorted(out)
def Wmat(a, T):
    L = 1 << a; S = range(-L, L); Cp0 = 1
    for q in S:
        if q: Cp0 *= -q
    CT = 1
    for q in S: CT *= (T - q)
    return [[(Fraction(CT, Cp0 * (T - (y - x)))) for y in range(L)] for x in range(L)]
def tm(fr): return (fr.numerator * pow(fr.denominator, -1, MOD)) % MOD
def bits(m, n): return [(m >> j) & 1 for j in range(n)]
def oneway(W, rows, cols):
    n = len(cols); out = []
    for Smask in range(1 << len(rows)):
        S = [rows[k] for k in range(len(rows)) if (Smask >> k) & 1]; row = []
        for Xm in range(1 << n):
            xb = bits(Xm, n); val = 1
            for x in S:
                val = val * sum(tm(W[x][c]) for k, c in enumerate(cols) if xb[k]) % MOD
            row.append(val)
        out.append(row)
    return out
def flattening(W, L, rg, cg, rh, ch):
    h = L // 2; Wm = [[tm(e) for e in r] for r in W]; A = []
    idx = list(itertools.product(range(1 << h), repeat=2))
    for (Yg, Xg) in idx:
        yg = bits(Yg, h); xg = bits(Xg, h); row = []
        for (Yh, Xh) in idx:
            yh = bits(Yh, h); xh = bits(Xh, h); X = [0] * L; Y = [0] * L
            for k, c in enumerate(cg): X[c] = xg[k]
            for k, c in enumerate(ch): X[c] = xh[k]
            for k, r in enumerate(rg): Y[r] = yg[k]
            for k, r in enumerate(rh): Y[r] = yh[k]
            val = 1
            for x in range(L):
                s = sum(Wm[x][y] for y in range(L) if X[y]) % MOD
                val = val * (s if Y[x] else (1 - s)) % MOD
                if val == 0: break
            row.append(val)
        A.append(row)
    return A
def cnt(vals, p): return sum(1 for v in vals if v < p)
PS = (2, 4, 6, 8, 10, 12, 16, 20, 24, 32, 40)
for (a, i) in [(2, 0), (2, 1), (3, 0), (3, 1), (3, 2)]:
    L = 1 << a; t = 1 << i; T = t + (1 << (a + 1)); W = Wmat(a, T)
    rg = [x for x in range(L) if (x >> i) & 1 == 0]; cg = [y for y in range(L) if (y >> i) & 1 == 1]
    rh = [x for x in range(L) if (x >> i) & 1 == 1]; ch = [y for y in range(L) if (y >> i) & 1 == 0]
    M1 = smith_vals(oneway(W, rg, ch), 60); M2 = smith_vals(oneway(W, rh, cg), 60)
    tens = sorted(u + v for u in M1 for v in M2)
    F = smith_vals(flattening(W, L, rg, cg, rh, ch), 41)
    print(f"a={a} i={i}: M1={M1}\n   M2={M2}")
    print("   R_p(F) / R_p(M1 x M2): " + " ".join(f"p={p}:{cnt(F,p)}/{cnt(tens,p)}" for p in PS)); sys.stdout.flush()
