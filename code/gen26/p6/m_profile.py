# Smith profile of the one-direction Boolean product matrix M[S,X'] = prod_{x in S} (E X')_x, E = crossing block.
# Compare with: tensor over tiles of M_{a-i,0}; full Smith-monomial weights; square-free (subset-sum) weights.
import sys, itertools
sys.setrecursionlimit(10000)
from fractions import Fraction
exec(open(sys.argv[1]).read().split("def cnt")[0].split("PS =")[0]) if False else None
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
                    if best is None or v < best[0]: best = (v, r, c)
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
def Eblock(a, t, zeta):
    L = 1 << a; S = range(-L, L); T = t + (1 << (a + 1)) * zeta; Cp0 = 1
    for q in S:
        if q: Cp0 *= -q
    CT = 1
    for q in S: CT *= (T - q)
    i = v2(t); B = [x for x in range(L) if (x >> i) & 1 == 0]
    f = lambda fr: (fr.numerator * pow(fr.denominator, -1, MOD)) % MOD
    return [[f(Fraction(CT, Cp0 * (T - (y - x)))) for y in B] for x in B]
def Mmat(E):
    n = len(E); out = []
    lin = [[0] * (1 << n) for _ in range(n)]
    for x in range(n):
        for Xm in range(1, 1 << n):
            low = Xm & -Xm; k = low.bit_length() - 1
            lin[x][Xm] = (lin[x][Xm ^ low] + E[x][k]) % MOD
    for Smask in range(1 << n):
        row = [1] * (1 << n)
        for x in range(n):
            if (Smask >> x) & 1:
                lx = lin[x]; row = [r * l % MOD for r, l in zip(row, lx)]
        out.append(row)
    return out
def sform(a, i): return [(a - i + 1) + 2 * k + 2 * sum(v2(m) for m in range(1, k + 1)) for k in range(1 << (a - i - 1))] * (1 << i)
def monos(ws, cap, sqfree):
    res = [0]
    for w in ws:
        new = []
        for r in res:
            e = 1
            while r + e * w < cap:
                new.append(r + e * w); e += 1
                if sqfree: break
        res = res + new
    return sorted(res)
prof = {}
cases = [(1, 0), (2, 0), (3, 0), (4, 0), (3, 1), (4, 1), (4, 2), (4, 3)]
for (a, i) in cases:
    M = Mmat(Eblock(a, 1 << i, 1)); cap = 400
    s = smith_vals(M, cap); prof[(a, i)] = s
    msg = f"a={a} i={i} n={len(s)} Smith(M)={s}"
    if i > 0 and (a - i, 0) in prof:
        tens = [0]
        for _ in range(1 << i): tens = sorted(u + v for u in tens for v in prof[(a - i, 0)])
        msg += f"\n   tile-tensor conj: {tens == s}"
    ws = sform(a, i); full = monos(ws, max(s) + 1, False)[:len(s)]; sq = monos(ws, 10**6, True)
    le_sq = all(u <= v for u, v in zip(s, sq)); ge_full = all(u >= v for u, v in zip(s, full))
    msg += f"\n   square-free sorted={sq}\n   termwise Smith<=sqfree: {le_sq};  Smith>=full-monomial: {ge_full}"
    print(msg); sys.stdout.flush()
