# Exact "true pairing rank" evaluator for lattice vertex algebras V_L and modules V_{lam+L}, in a regraded grading
#   g_xi = L0 - xi(0):  weight of e^p (x) (oscillators of level k) is |p|^2/2 - <xi,p> + k.
# Objects (all exact, integer / rational arithmetic, ranks with python-flint):
#   B_g(a,b)  = coefficient of e^g in a(H-1)b   (H = wt0 a + wt0 b - |g|^2/2), computed by Wick contraction
#               (3-point function <e^g| Y(a,1) b>), cocycle signs dropped (they rescale whole blocks: ranks unchanged).
#   O_g       = o_xi(e^g) : M^{xi}_{E0} -> M^{xi}_{E0}, the g_xi-weight-preserving mode of e^g, computed from
#               Y(e^g,z) = exp(sum g_{-n} z^n/n) exp(-sum g_n z^-n/n) e^g z^{g(0)}.
#   R_g(ha,hb) = rank of B_g on V^{xi}_{ha} x V^{xi}_{hb};   Lab_g(E0) = rank O_g on M^{xi}_{E0}.
# A construction (A1)-(A4) with K = 0 (c_xv in V^xi_0 = span e^g, g in S_xi) satisfies
#   d [x=v][y=u] = sum_g B_g(a_x,b_v) <u|O_g|y>   =>   N^2 <= sum_g rk(B_g|) rk(O_g|)  and  N <= sum_g min(R_g, Lab_g).
import itertools, math, sys
from fractions import Fraction
from functools import lru_cache
import flint
import numpy as np

def ip(G, u, v):  # <u,v> with u,v in lattice-basis coordinates (rationals ok)
    return sum(u[i]*G[i][j]*v[j] for i in range(len(u)) for j in range(len(v)))

@lru_cache(maxsize=None)
def fock(r, k):
    """basis monomials of level k for rank r: sorted tuples of (c,m), m>=1, sum m = k"""
    out = []
    def rec(rem, maxpart, cur):
        if rem == 0: out.append(tuple(cur)); return
        for m in range(min(rem, maxpart[0]), 0, -1):
            for c in range(r):
                if (m, c) > maxpart: continue
                rec(rem - m, (m, c), cur + [(m, c)])
    rec(k, (k, r), [])
    return [tuple(sorted(((c, m) for (m, c) in t))) for t in out]

def lattice_points(G, R2max, coset=None):
    """points p = coset + n (n integer vector) with |p|^2 <= R2max (coset in basis coords, rational)"""
    r = len(G)
    coset = coset or [0]*r
    Gm = np.array(G, float); ev = np.linalg.eigvalsh(Gm).min()
    B = int(math.ceil(math.sqrt(R2max/ev))) + 2
    pts = []
    for n in itertools.product(range(-B, B+1), repeat=r):
        p = tuple(Fraction(coset[i]) + n[i] for i in range(r))
        if ip(G, p, p) <= R2max: pts.append(p)
    return pts

def wt(G, xi, p, k=0):
    return Fraction(ip(G, p, p), 2) - ip(G, xi, p) + k

# ---------- pairing form B_g ----------
def Bval(G, beta, amono, delta, bmono):
    """<e^{beta+delta}| Y(a,1) b>, a = e^beta (x) prod alpha^{c}_{-m}, b = e^delta (x) bmono (cocycle dropped)"""
    r = len(G)
    Gd = [sum(G[c][j]*delta[j] for j in range(r)) for c in range(r)]
    Gb = [sum(G[c][j]*beta[j] for j in range(r)) for c in range(r)]
    poly = {tuple(bmono): Fraction(1)}
    for (c, m) in amono:
        new = {}
        sgn = (-1)**(m-1)
        for mono, co in poly.items():
            if Gd[c]:
                new[mono] = new.get(mono, 0) + sgn*Gd[c]*co
            done = set()
            for idx, (d, l) in enumerate(mono):
                if (d, l) in done or G[c][d] == 0: continue
                done.add((d, l))
                mult = mono.count((d, l))
                lst = list(mono); lst.remove((d, l)); key = tuple(lst)
                new[key] = new.get(key, 0) + sgn*math.comb(l+m-1, m-1)*l*G[c][d]*mult*co
        poly = {k2: v for k2, v in new.items() if v}
    tot = Fraction(0)
    for mono, co in poly.items():
        t = co
        for (d, l) in mono: t *= -Gb[d]
        tot += t
    return tot

def rank_q(rows):
    if not rows or not rows[0]: return 0
    return flint.fmpq_mat([[flint.fmpq(x.numerator, x.denominator) for x in row] for row in rows]).rank()

def R_channel(G, xi, g, ha, hb, verbose=False):
    """rank of B_g on V^xi_ha x V^xi_hb (sum over blocks beta, delta=g-beta)"""
    r = len(G)
    tot = 0; blocks = []
    R2 = 2*(max(ha, hb) + 2*ip(G, xi, xi) + 4)
    for beta in lattice_points(G, R2):
        delta = tuple(g[i]-beta[i] for i in range(r))
        ka = ha - wt(G, xi, beta); kb = hb - wt(G, xi, delta)
        if ka < 0 or kb < 0 or ka.denominator != 1 or kb.denominator != 1: continue
        A = fock(r, int(ka)); Bs = fock(r, int(kb))
        M = [[Bval(G, beta, am, delta, bm) for bm in Bs] for am in A]
        rk = rank_q(M); tot += rk
        if rk: blocks.append((beta, int(ka), int(kb), len(A), len(Bs), rk))
    if verbose: print("   blocks (beta, ka, kb, dimA, dimB, rank):", blocks)
    return tot, blocks

# ---------- label map O_g ----------
def shift_poly(mono, s):  # substitute x_{c,n} -> x_{c,n} - s_c
    poly = {(): Fraction(1)}
    for (c, n) in mono:
        new = {}
        for m, co in poly.items():
            k1 = tuple(sorted(m + ((c, n),))); new[k1] = new.get(k1, 0) + co
            if s[c]:
                new[m] = new.get(m, 0) - s[c]*co
        poly = {k: v for k, v in new.items() if v}
    return poly

@lru_cache(maxsize=None)
def exp_part(g, r, j):
    """level-j part of exp(sum_n Gamma_n/n), Gamma_n = sum_c g_c x_{c,n}; as dict mono->coef"""
    out = {}
    def parts(n, maxp):
        if n == 0: yield (); return
        for k in range(min(n, maxp), 0, -1):
            for rest in parts(n-k, k): yield (k,) + rest
    for lam in parts(j, j):
        poly = {(): Fraction(1)}
        from collections import Counter
        for n, mult in Counter(lam).items():
            fac = Fraction(1, n**mult * math.factorial(mult))
            for _ in range(mult):
                new = {}
                for m, co in poly.items():
                    for c in range(r):
                        if g[c]:
                            key = tuple(sorted(m + ((c, n),))); new[key] = new.get(key, 0) + co*g[c]
                poly = new
            poly = {k: v*fac for k, v in poly.items()}
        for k, v in poly.items(): out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v}

def O_apply(G, g, p, mono, kout):
    r = len(G)
    Gg = [sum(G[c][j]*g[j] for j in range(r)) for c in range(r)]
    sp = shift_poly(mono, Gg)
    out = {}
    for m, co in sp.items():
        j = kout - sum(n for (_, n) in m)
        if j < 0: continue
        for m2, c2 in exp_part(tuple(g), r, j).items():
            key = tuple(sorted(m + m2)); out[key] = out.get(key, 0) + co*c2
    return {k: v for k, v in out.items() if v}

def Lab_channel(G, xi, g, E0, coset=None):
    r = len(G); tot = 0; blocks = []
    for p in lattice_points(G, 2*(E0 + 2*ip(G, xi, xi) + 4), coset):
        k = E0 - wt(G, xi, p)
        q = tuple(p[i]+g[i] for i in range(r)); k2 = E0 - wt(G, xi, q)
        if k < 0 or k2 < 0 or k.denominator != 1: continue
        dom = fock(r, int(k)); cod = fock(r, int(k2)); idx = {m: i for i, m in enumerate(cod)}
        rows = []
        for mono in dom:
            v = [Fraction(0)]*len(cod)
            for m2, c2 in O_apply(G, g, p, mono, int(k2)).items(): v[idx[m2]] += c2
            rows.append(v)
        rk = rank_q(rows); tot += rk; blocks.append((p, int(k), int(k2), len(dom), len(cod), rk))
    return tot, blocks

def dimM(G, xi, E, coset=None):
    r = len(G); tot = 0
    for p in lattice_points(G, 2*(E + 2*ip(G, xi, xi) + 4), coset):
        k = E - wt(G, xi, p)
        if k >= 0 and k.denominator == 1: tot += len(fock(r, int(k)))
    return tot

def S_xi(G, xi):
    return [g for g in lattice_points(G, 4*ip(G, xi, xi) + 1) if wt(G, xi, g) == 0]

if __name__ == "__main__":
    pass
