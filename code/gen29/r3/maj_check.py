# G29-R3 independent check of the Majorana fact and family (own code, not reusing H1 scripts).
# F2^{2k}: w = c | (d << k).  Q(w) = c.d, B(w,v) = c_w.d_v + d_w.c_v.  beta = (-1)^Q.
# Frames: RANDOM 2k pairwise-B=1 vectors (backtracking), completed by their sum (not Jordan-Wigner).
import sys, random
from fractions import Fraction
from math import comb
import numpy as np
random.seed(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
k = int(sys.argv[1]); n2 = 4**k; msk = (1 << k) - 1
pc = np.vectorize(lambda x: bin(int(x)).count('1'))
W = np.arange(n2, dtype=np.int64)
POP = np.array([bin(i).count('1') for i in range(n2)], dtype=np.int64)
def Q(x): return POP[(x & msk) & (x >> k)] & 1
def B(x, v): return (POP[(x & msk) & (v >> k)] + POP[(x >> k) & (v & msk)]) & 1
def frame():
    while True:
        S = []
        cand = list(range(1, n2)); random.shuffle(cand)
        def bt(S, start):
            if len(S) == 2*k: return S
            for idx in range(start, len(cand)):
                v = cand[idx]
                if all(B(np.int64(v), np.int64(s)) == 1 for s in S):
                    r = bt(S+[v], idx+1)
                    if r: return r
                if len(S) == 0 and idx > 50: return None
            return None
        r = bt([], 0)
        if r: break
    t = 0
    for v in r: t ^= v
    return r + [t]
fr = frame()
assert all(B(np.int64(a), np.int64(b)) == 1 for i, a in enumerate(fr) for b in fr[i+1:]), 'not pairwise 1'
assert np.bitwise_xor.reduce(np.array(fr)) == 0
beta = 1 - 2*Q(W)
Ylin = np.stack([B(W, np.int64(s)) for s in fr])          # linear forms
Yaff = np.stack([(Q(np.int64(s)) + B(W, np.int64(s))) & 1 for s in fr])  # affine forms y_s = Q(s)+B(s,w)
def isfunc(J):
    ok = True; F = {}
    for j in range(4):
        m = (J % 4 == j)
        if m.any():
            vals = set(beta[m].tolist())
            if len(vals) > 1: ok = False
            else: F[j] = vals.pop()
    return ok, F
Jl = Ylin.sum(0); Ja = Yaff.sum(0)
okl, _ = isfunc(Jl); oka, F = isfunc(Ja)
# injectivity onto parity class
codes = (Yaff * (1 << np.arange(2*k+1))[:, None]).sum(0)
inj = len(set(codes.tolist())) == n2
par = int(Ja[0] % 2); assert np.all(Ja % 2 == par)
print(f'k={k} frame Q-values={[int(Q(np.int64(s))) for s in fr]} | linear forms: beta=F(|y| mod4)? {okl} | affine: {oka} F={F} par={par} injective={inj}', flush=True)
if len(sys.argv) > 3 and sys.argv[3] == 'fam':
    u = (Ja - par)//2; umax = (2*k+1-par)//2
    Fu = lambda x: F[(2*x+par) % 4]
    # symplectic Fourier: c(s) = N^{-2} sum_v phi(v) (-1)^{B(s,v)}; B(s,v) = <s, Jv> with J swapping halves -> plain WHT of phi after relabel
    def wht(a):
        a = a.copy(); h = 1
        while h < len(a):
            a = a.reshape(-1, 2*h); x = a[:, :h].copy(); y = a[:, h:].copy(); a[:, :h] = x+y; a[:, h:] = x-y; a = a.reshape(-1); h *= 2
        return a
    res = []
    for lo in range(umax+1):
        for hi in range(lo, umax+1):
            win = range(lo, hi+1)
            def f(x):
                s = Fraction(0)
                for i in win:
                    t = Fraction(Fu(i))
                    for j in win:
                        if j != i: t *= Fraction(x-j, i-j)
                    s += t
                return s
            fv = [f(x) for x in range(umax+1)]
            den = 1
            for v in fv: den = den*v.denominator//np.gcd(den, v.denominator)
            phi_int = np.array([int(v*den) for v in fv], dtype=object)[u]
            hat = wht(np.array(phi_int, dtype=object))
            q = int(sum(1 for h in hat if h != 0))
            p = int(sum(1 for i in range(len(u)) if fv[u[i]] != Fu(u[i])))
            d = hi - lo
            qf = sum(comb(2*k+1, j) for j in range(d+1)) if d <= k else None
            pf = sum(comb(2*k+1, 2*x+par) for x in range(umax+1) if x < lo or x > hi)
            res.append((q, p, lo, hi, qf, pf))
    bad = [r for r in res if (r[4] is not None and r[0] > r[4]) or r[1] != r[5]]
    eqq = sum(1 for r in res if r[4] is not None and r[0] == r[4])
    print(f'  windows={len(res)} q<=formula & p==formula violations: {bad[:5]} (#{len(bad)}); q==formula in {eqq}/{sum(1 for r in res if r[4] is not None)}')
    res.sort(); fr_ = []; bp = 10**30
    for r in res:
        if r[1] < bp: fr_.append(r[:4]); bp = r[1]
    print('  frontier (q,p,lo,hi):', fr_)
    print('  balanced min max(p,q):', min(max(r[0], r[1]) for r in res), flush=True)
