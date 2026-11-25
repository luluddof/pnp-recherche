# Majorana Bell construction (H1): frame of 2k+1 vectors pairwise B=1, f = F(J); window interpolation. Returns g (function) and support.
import numpy as np, itertools
from core import *
def Bform(k):
    # f(z) = (-1)^{x.y}; polar form B(z,w) = x.y' + y.x'
    def B(z, w):
        m = (1 << k) - 1
        return (popc((z & m) & (w >> k)) + popc((z >> k) & (w & m))) & 1
    return B
def frame(k, rng=None):
    # greedy/backtracking search for 2k+1 vectors pairwise B=1, spanning
    B = Bform(k); D = 4**k
    def rec(cur):
        if len(cur) == 2*k + 1: return cur
        start = cur[-1] + 1 if cur else 1
        for v in range(start, D):
            if all(B(v, u) == 1 for u in cur):
                r = rec(cur + [v])
                if r: return r
        return None
    return rec([])
def majorana_g(k, d, w0):
    D = 4**k; f = beta(k); B = Bform(k); S = frame(k); z = np.arange(D)
    Y = np.array([[B(s, int(w)) for s in S] for w in z])      # bits B(s_i,w); Q(s)-shift handled by fitting
    # f is a function of y (affine forms): fit via c: f(w) = F(|y(w)+c|); find c in {0,1}^{2k+1} making f a function of weight
    for c in itertools.product([0, 1], repeat=2*k+1):
        J = ((Y + np.array(c)) % 2).sum(1); ok = True; F = {}
        for j, val in zip(J, f):
            if F.setdefault(j, val) != val: ok = False; break
        if ok: break
    Js = sorted(F); par = Js[0] % 2
    u = (J - par) // 2; us = sorted(set(u))
    W = [w0 + i for i in range(d + 1)]
    # degree-d polynomial in u interpolating F on window W
    coef = np.polyfit(W, [F[2*x + par] for x in W], d)
    g = np.polyval(coef, u)
    gh = fwht(g) / D; supp = np.abs(gh) > 1e-9
    p = int((np.abs(g - f) > 1e-9).sum())
    return g, supp, p, J
if __name__ == '__main__':
    for k in (2, 3):
        for d in (0, 1, 2):
            for w0 in range(0, k + 1 - d):
                g, supp, p, J = majorana_g(k, d, w0)
                print('k', k, 'd', d, 'w0', w0, 'q', supp.sum(), 'p', p)
