"""G18c: 2-adic rank profile of OpenAI's kernel on one packet.
For node t and an evaluation T = t + alpha*zeta (alpha = 2^{a+1}), the flattened kernel K(T) in direction
i = v2(t) is a matrix over Z_(2) (Lemma D.4). Over Z/2^p the minimal number of rank-one terms of a matrix
is #{Smith invariants with v2 < p}.  We print N_p = that number for p = 1, 2, ... : the exact
'DF-rank mod 2^p' of the one-packet kernel, to compare with 2^L (rank over Q) and with Claude's
weight-window profile dim W_{<=x}."""
import sys, itertools
from fractions import Fraction
import flint
a = int(sys.argv[1]); L = 1 << a; alpha = 1 << (a + 1)
zetas = [int(z) for z in sys.argv[2].split(',')] if len(sys.argv) > 2 else [1]
S = list(range(-L, L))
Cp0 = 1
for t in S:
    if t: Cp0 *= -t
def Wmat(T):
    T = Fraction(T)
    M = [[None] * L for _ in range(L)]
    for x in range(L):
        for y in range(L):
            v = Fraction(1, Cp0)
            for t in S:
                if t != y - x: v *= (T - t)
            M[x][y] = v
    return M
bits = list(itertools.product([0, 1], repeat=L))
def v2int(n):
    if n == 0: return 10**9
    n = abs(n); k = 0
    while n % 2 == 0: n //= 2; k += 1
    return k
def v2(t): return v2int(t)
def kernel(T):
    Wm = Wmat(T)
    K = {}
    for xi, X in enumerate(bits):
        forms = [sum((Wm[x][y] for y in range(L) if X[y]), Fraction(0)) for x in range(L)]
        for yi, Y in enumerate(bits):
            f = Fraction(1)
            for x in range(L):
                f *= forms[x] if Y[x] else (1 - forms[x])
                if f == 0: break
            K[(yi, xi)] = f
    return K
def flatten(K, i):
    g0 = [('Y', x) for x in range(L) if not (x >> i) & 1] + [('X', y) for y in range(L) if (y >> i) & 1]
    g1 = [('Y', x) for x in range(L) if (x >> i) & 1] + [('X', y) for y in range(L) if not (y >> i) & 1]
    M = [[Fraction(0)] * (1 << len(g1)) for _ in range(1 << len(g0))]
    for yi, Y in enumerate(bits):
        for xi, X in enumerate(bits):
            v = K[(yi, xi)]
            if v == 0: continue
            r = sum(((Y[p] if s == 'Y' else X[p]) << k) for k, (s, p) in enumerate(g0))
            c = sum(((Y[p] if s == 'Y' else X[p]) << k) for k, (s, p) in enumerate(g1))
            M[r][c] = v
    return M
for t in S:
    if t == 0 or v2(t) >= a: continue
    i = v2(t)
    for zeta in zetas:
        K = kernel(t + alpha * zeta)
        M = flatten(K, i)
        den = 1
        for row in M:
            for v in row:
                d = v.denominator
                den = den * d // __import__('math').gcd(den, d)
        assert den % 2 == 1, "entry not 2-integral (Lemma D.4 violated?)"
        Mi = flint.fmpz_mat([[int(v * den) for v in row] for row in M])
        snf = Mi.snf()
        inv = sorted(v2int(int(snf[k, k])) for k in range(min(snf.nrows(), snf.ncols())))
        hist = {}
        for e in inv: hist[e] = hist.get(e, 0) + 1
        cum = []; c = 0
        for e in sorted(hist):
            c += hist[e]; cum.append((e, c))
        print(f"a={a} L={L} t={t:+d} dir={i} zeta={zeta}: rank_Q={sum(1 for e in inv if e<10**9)} of {len(inv)};"
              f" Smith 2-adic valuations: histogram {dict(sorted(hist.items()))}; N_p (rank mod 2^p) = "
              + ", ".join(f"p>{e}:{c}" for e, c in cum if e < 10**9), flush=True)
