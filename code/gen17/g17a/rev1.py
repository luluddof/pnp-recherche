#!/usr/bin/env python3
"""G17a independent checks (reviewer code, not derived from the prover's checks.py)."""
from fractions import Fraction as Fr
import itertools, math, random
INF = 10**9
def v2(x):
    x = Fr(x)
    if x == 0: return INF
    n, d, k = abs(x.numerator), x.denominator, 0
    while n % 2 == 0: n //= 2; k += 1
    while d % 2 == 0: d //= 2; k -= 1
    return k
# ---- polynomials as coefficient lists (low degree first)
def pmul(A, B):
    R = [Fr(0)] * (len(A) + len(B) - 1)
    for i, x in enumerate(A):
        if x == 0: continue
        for j, y in enumerate(B): R[i + j] += x * y
    return R
def padd(A, B):
    n = max(len(A), len(B)); return [(A[i] if i < len(A) else 0) + (B[i] if i < len(B) else 0) for i in range(n)]
def pscale(A, c): return [c * x for x in A]
def ppow(A, k):
    R = [Fr(1)]
    for _ in range(k): R = pmul(R, A)
    return R
def peval(A, x):
    r = Fr(0)
    for c in reversed(A): r = r * x + c
    return r
def pcompose_affine(A, t, al):
    """A(t + al*z) as polynomial in z"""
    R = [Fr(0)]; lin = [Fr(t), Fr(al)]
    for c in reversed(A): R = padd(pmul(R, lin), [c])
    return R
def gauss(A): return min(v2(c) for c in A)
def deg(A):
    d = len(A) - 1
    while d > 0 and A[d] == 0: d -= 1
    return d

def grid(a):
    L = 2 ** a; return L, list(range(-L, L))

# ---------- R1: the Cauchy deformation -------------------------------------
def R1(a):
    L, S = grid(a); al = 2 ** (a + 1)
    C = [Fr(1)]
    for t in S: C = pmul(C, [Fr(-t), Fr(1)])
    def Cprime(t):
        p = Fr(1)
        for q in S:
            if q != t: p *= (t - q)
        return p
    ok = all(v2(Cprime(t)) == v2(Cprime(0)) for t in S)
    def Wpoly(x, y):
        p = [Fr(1)]
        for q in S:
            if q != y - x: p = pmul(p, [Fr(-q), Fr(1)])
        return pscale(p, 1 / Cprime(0))
    Wp = {(x, y): Wpoly(x, y) for x in range(L) for y in range(L)}
    ok &= all(deg(p) == 2 * L - 1 for p in Wp.values())
    for (x, y), p in Wp.items():
        ok &= peval(p, 0) == (1 if x == y else 0)
        for t in S:
            if t == 0: continue
            w = peval(p, t)
            ok &= (w != 0) == (y - x == t)
            if w != 0:
                ok &= (w == Cprime(t) / Cprime(0)) and v2(w) == 0
                i = v2(t)
                ok &= ((x >> i) & 1) != ((y >> i) & 1)      # NOTE: lower bits need not agree in general? test:
            # Lemma D.4: W(t+al z) integral; entries with y-x != t divisible by 2 (indeed 2^{a+1-v2(t-(y-x))})
            g = gauss(pcompose_affine(p, t, al))
            if y - x == t: ok &= (g == 0)
            else: ok &= (g == a + 1 - v2(t - (y - x)))
            # reviewer's extra: integrality already on the larger disc t + 2^a z'
            ok &= gauss(pcompose_affine(p, t, 2 ** a)) >= 0
    return ok
print("R1  W: deg 2L-1, W(0)=I, W(t)=unit*shift_t, bit v2(t) flips, Lemma D.4 valuations, integrality on disc radius 2^-a:",
      [R1(a) for a in (1, 2, 3)])

# ---------- R2: identity (D.10), arbitrary profiles (also NON valuation-dependent: brute definition only)
def loss_direct(S, q, t):
    return -sum(q[s] * (v2(s) - v2(t - s)) for s in S if s not in (0, t))
def lam_formula(a, prof, i):
    return sum(Fr(j, 2 ** j) * (prof[i] - prof[i + j]) for j in range(1, a - i + 1))
random.seed(7)
ok = True
for a in range(1, 8):
    L, S = grid(a)
    for _ in range(5):
        prof = [random.randint(0, 12) for _ in range(a + 1)]
        q = {t: prof[v2(t)] for t in S if t != 0}
        for t in S:
            if t == 0: continue
            i = v2(t)
            ok &= Fr(loss_direct(S, q, t)) == Fr(L, 2 ** i) * lam_formula(a, prof, i)
print("R2  (D.10): -v2(E_t(0)) = D_t * sum_j j/2^j (q(i)-q(i+j)) for random valuation profiles:", ok)

# ---------- R3: Theorem B inequality, exhaustive/random + LP
def maxlam(a, u, prof):
    used = [i for i in range(u) if prof[i] > 0]
    return max(lam_formula(a, prof, i) for i in used)
bad = 0; worst = None
for a in range(1, 9):
    L = 2 ** a
    for u in range(1, a + 1):
        for _ in range(3000):
            prof = [random.choice([0, 0, 1, 2, 3, 5, 8, 13, 40]) for _ in range(u)] + [0] * (a + 2)
            tot = sum(Fr(L, 2 ** i) * prof[i] for i in range(u))   # = sum_t q_t
            if tot == 0: continue
            # normalise: M := largest real with tot > M(2L-1)  -> M = tot/(2L-1)
            Mreal = tot / (2 * L - 1)
            ratio = maxlam(a, u, prof) / ((2 - Fr(1, L)) * Mreal / (4 * u))
            if ratio < 1: bad += 1
            if worst is None or ratio < worst[0]: worst = (ratio, a, u, prof[:u])
print("R3  Thm B inequality max_i Lambda_i >= (2-1/L)M/(4u): violations =", bad, " smallest ratio seen:", float(worst[0]), worst[1:])
try:
    import numpy as np
    from scipy.optimize import linprog
    print("R3b LP: min over profiles of u*max_i Lambda_i (M=1), and paper taper (h=3/u), minimal taper (h=1/u):")
    for a, u in [(4, 2), (8, 4), (16, 8), (24, 12), (40, 20), (60, 30), (30, 30)]:
        L = 2.0 ** a
        c = np.zeros(u + 1); c[-1] = 1; A = []; b = []
        for i in range(u):
            row = np.zeros(u + 1)
            for j in range(1, a - i + 1):
                row[i] += j / 2.0 ** j
                if i + j < u: row[i + j] -= j / 2.0 ** j
            row[-1] = -1; A.append(row); b.append(0)
        row = np.zeros(u + 1)
        for i in range(u): row[i] = -2.0 ** (-i)
        A.append(row); b.append(-(2 - 1 / L))
        r = linprog(c, A_ub=np.array(A), b_ub=np.array(b), bounds=[(0, None)] * (u + 1))
        def tl(h):
            qq = [h * (u - i) for i in range(u)] + [0] * (a + 2)
            return max(sum(j / 2.0 ** j * (qq[i] - qq[i + j]) for j in range(1, a - i + 1)) for i in range(u))
        print(f"      a={a:3d} u={u:3d}  LP*u={r.fun*u:7.4f}  paper*u={tl(3.0/u)*u:7.3f}  minimal-taper*u={tl(1.0/u)*u:7.3f}  proved*u={(2-1/L)/4:.4f}")
except Exception as e:
    print("R3b skipped", e)

# ---------- R4: tuple counts (own recurrence), weights w>=1, n variables per weight, sum w*beta < p
def count(n, p, wmin=1, wmax=None):
    coef = [0] * p; coef[0] = 1
    for w in range(wmin, p if wmax is None else min(p, wmax + 1)):
        for _ in range(n):
            for k in range(w, p): coef[k] += coef[k - w]
    return sum(coef)
print("R4  log2 #tuples vs prover floor sqrt(np)-n ; also with weights shifted by s (true divisibility b+1+a-i) and capped at D/2:")
for n, p, s, cap in [(2, 100, 0, None), (4, 400, 0, None), (4, 400, 6, None), (4, 400, 0, 40), (8, 1600, 8, 64), (16, 256, 0, None), (40, 30, 0, None)]:
    c = count(n, p, 1 + s, None if cap is None else cap + s)
    print(f"      n={n:3d} p={p:5d} shift={s} cap={cap}:  log2={math.log2(c):8.2f}  floor={math.sqrt(n*p)-n:8.2f}  Meinardus={2*math.sqrt(math.pi**2/6*n*p)/math.log(2):8.2f}")
