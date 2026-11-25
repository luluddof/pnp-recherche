# Majorana Bell windows (after G29-H1): phi on w = x | (z<<k), approximating beta(w) = (-1)^{x.z}.
# returns list of dicts: q (Fourier sparsity), p (#disagreements), phiL (integer array = L*phi), L, supp (WHT support)
import numpy as np
from fractions import Fraction
from math import lcm
from core import popc, fwht
def majorana_windows(k, frontier_only=True):
    n = 4**k; msk = (1 << k) - 1
    A = []
    for j in range(k):
        for P in 'XY':
            a = 1 << j; b = ((1 << j) - 1) | ((1 << j) if P == 'Y' else 0); A.append(a | (b << k))
    t = 0
    for v in A: t ^= v
    A.append(t)
    w = np.arange(n)
    Qv = lambda x: popc((x & msk) & (x >> k)) & 1
    Bv = lambda x, s: (popc((x & msk) & (s >> k)) + popc((s & msk) & (x >> k))) & 1
    J = np.zeros(n, np.int64)
    for s in A: J += (Qv(np.array([s]))[0] + Bv(w, s)) & 1
    beta = 1 - 2*Qv(w); par = int(J[0] % 2)
    F = {}
    for j4 in (par, par+2):
        vals = set(beta[J % 4 == j4]); assert len(vals) == 1; F[j4] = int(vals.pop())
    u = (J - par)//2; Fu = lambda uu: F[(2*uu + par) % 4]
    out = []
    for lo in range(k+1):
        for hi in range(lo, k+1):
            W = list(range(lo, hi+1))
            def f(x):
                s = Fraction(0)
                for i in W:
                    tt = Fraction(Fu(i))
                    for jj in W:
                        if jj != i: tt *= Fraction(x - jj, i - jj)
                    s += tt
                return s
            fv = [f(x) for x in range(k+1)]
            L = 1
            for v in fv: L = lcm(L, v.denominator)
            fi = np.array([int(v*L) for v in fv], dtype=object)
            phiL = fi[u]
            p = int(sum(np.sum(u == x) for x in range(k+1) if fv[x] != Fu(x)))
            wh = fwht(phiL); supp = [i for i in range(n) if wh[i] != 0]
            out.append(dict(q=len(supp), p=p, phiL=phiL, L=L, supp=supp, win=(lo, hi), wh=wh))
    # add exact beta (q = N^2 ... ) and zero
    if frontier_only:
        out.sort(key=lambda d: (d['q'], d['p'])); fr = []; bp = 10**18
        for d in out:
            if d['p'] < bp: fr.append(d); bp = d['p']
        return fr, beta
    return out, beta
if __name__ == '__main__':
    import sys
    k = int(sys.argv[1]); fr, beta = majorana_windows(k)
    print([(d['q'], d['p'], d['win']) for d in fr])
