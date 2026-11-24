#!/usr/bin/env python3
"""G17a part 2: (a) Appendix D's kernel identity is of DF form on a toy packet (exact, over Q);
(b) Theorem A sanity on that toy; (c) product witnesses for the necessity of the precision;
(d) TRUE minimal precision (local Smith form over Z_(2)) for the black-box class and for the
    structured class (Z_(2)-span of the kernel entries K_T(Y,X))."""
from fractions import Fraction as Fr
import itertools, math, sys
from rev1_lib import *

def setup(a):
    L, S = grid(a); al = 2 ** (a + 1)
    def Cprime(t):
        p = Fr(1)
        for q in S:
            if q != t: p *= (t - q)
        return p
    def Wpoly(d):
        p = [Fr(1)]
        for q in S:
            if q != d: p = pmul(p, [Fr(-q), Fr(1)])
        return pscale(p, 1 / Cprime(0))
    Wd = {d: Wpoly(d) for d in range(-L + 1, L)}
    return L, S, al, Wd

def kernel(a, slots):
    """K_T(Y,X) for one packet whose real slots are `slots` (each slot = one wire)."""
    L, S, al, Wd = setup(a)
    M = len(slots); K = {}
    for X in itertools.product((0, 1), repeat=M):
        forms = []
        for x in slots:
            f = [Fr(0)]
            for k, y in enumerate(slots):
                if X[k]: f = padd(f, Wd[y - x])
            forms.append(f)
        for Y in itertools.product((0, 1), repeat=M):
            p = [Fr(1)]
            for k in range(M):
                p = pmul(p, forms[k] if Y[k] else padd([Fr(1)], pscale(forms[k], -1)))
            K[(Y, X)] = p
    return K

def jets(P, t, al, q):
    c = pcompose_affine(P, t, al); c = c + [Fr(0)] * q
    return c[:q]

def series_inv(E, q):
    inv = [Fr(0)] * q; inv[0] = 1 / E[0]
    for k in range(1, q):
        inv[k] = -sum(E[j] * inv[k - j] for j in range(1, min(k, len(E) - 1) + 1)) / E[0]
    return inv

def hermite_weights(a, prof):
    """linear functional: F(0) = sum_{t,j} c[t][j] * [z^j]F(t+al z)   (formula (D.9))"""
    L, S = grid(a); al = 2 ** (a + 1); c = {}
    for t in S:
        if t == 0: continue
        q = prof[v2(t)]
        if q == 0: continue
        E = [Fr(1)]
        for s in S:
            if s in (0, t): continue
            E = pmul(E, ppow([Fr(-s, t - s), Fr(1, t - s)], prof[v2(s)]))
        E0 = peval(E, 0); Einv = series_inv(pcompose_affine(E, t, al), q)
        # [F/E]_{<q}(z0) = sum_k z0^k sum_{j<=k} f_j Einv_{k-j}
        z0 = Fr(-t, al)
        c[t] = [E0 * sum(z0 ** k * Einv[k - j] for k in range(j, q)) for j in range(q)]
    return c

def local_smith_minprec(rows, target, colscale=None):
    """rows: list of vectors (data coordinates) over Z_(2); target: list of scalars (value of F(0)).
    Lattice = Z_(2)-span of the pairs.  Returns min P such that  data = 0 mod 2^P  (after column
    scaling by 2^colscale[c])  implies target = 0 mod 2.  Also returns rank."""
    R = [list(map(Fr, r)) for r in rows]; n = len(R[0])
    if colscale: R = [[r[c] * Fr(2) ** colscale[c] for c in range(n)] for r in R]
    w = list(map(Fr, target)); m = len(R); rank = 0; need = 0
    usedrow = [False] * m; usedcol = [False] * n
    while True:
        best = None
        for i in range(m):
            if usedrow[i]: continue
            for c in range(n):
                if usedcol[c] or R[i][c] == 0: continue
                v = v2(R[i][c])
                if best is None or v < best[0]: best = (v, i, c)
        if best is None: break
        v, i, c = best; usedrow[i] = True; usedcol[c] = True; rank += 1
        piv = R[i][c]
        for k in range(m):                       # row operations (2-integral multipliers)
            if k == i or R[k][c] == 0: continue
            f = R[k][c] / piv; assert v2(f) >= 0
            R[k] = [R[k][j] - f * R[i][j] for j in range(n)]; w[k] -= f * w[i]
        for j in range(n):                       # column operations: clear pivot row
            if j != c: R[i][j] = Fr(0)
        if v2(w[i]) < 0: return None, rank       # target not even integral: class badly posed
        if v2(w[i]) == 0: need = max(need, v + 1)
        # (if v2(w[i]) = e > 0 nothing is required for the residue mod 2)
    for i in range(m):
        if not usedrow[i]: assert w[i] == 0, "target not determined by data"
    return need, rank

if __name__ == "__main__":
    # ---------------- (a),(b) toy packet a=1, slots 0,1 ----------------
    for a, slots, prof in [(1, [0, 1], [6, 0, 0]), (2, [0, 1, 2], [9, 5, 0, 0]), (2, [0, 1, 2, 3], [12, 6, 0, 0])]:
        L, S = grid(a); al = 2 ** (a + 1); M = len(slots); N = M * (2 * L - 1)
        K = kernel(a, slots)
        okI = all(peval(P, 0) == (1 if Y == X else 0) for (Y, X), P in K.items())
        okdeg = max(deg(P) for P in K.values()) <= N
        c = hermite_weights(a, prof)
        okH = True
        for (Y, X), P in K.items():
            val = sum(cj * fj for t in c for cj, fj in zip(c[t], jets(P, t, al, prof[v2(t)])))
            okH &= (val == (1 if Y == X else 0))
        okInt = all(gauss(pcompose_affine(P, t, al)) >= 0 for P in K.values() for t in S)
        # locality at nodes: with i=v2(t), group of row slot x = bit i of x, group of column slot y = 1-bit i of y
        okloc = True; termrank = 0
        import sympy
        for t in c:
            i = v2(t)
            rowbits = [k for k, x in enumerate(slots)]; 
            g0 = ([k for k, x in enumerate(slots) if (x >> i) & 1 == 0], [k for k, y in enumerate(slots) if (y >> i) & 1 == 1])
            g1 = ([k for k, x in enumerate(slots) if (x >> i) & 1 == 1], [k for k, y in enumerate(slots) if (y >> i) & 1 == 0])
            def idx(Y, X, g): return tuple(Y[k] for k in g[0]) + tuple(X[k] for k in g[1])
            r0 = sorted({idx(Y, X, g0) for (Y, X) in K}); r1 = sorted({idx(Y, X, g1) for (Y, X) in K})
            for j in range(prof[i]):
                Mx = sympy.zeros(len(r0), len(r1))
                for (Y, X), P in K.items():
                    Mx[r0.index(idx(Y, X, g0)), r1.index(idx(Y, X, g1))] = sympy.Rational(*(lambda f: (f.numerator, f.denominator))(jets(P, t, al, prof[i])[j]))
                rk = Mx.rank(); termrank += rk
                if j == 0: okloc &= (rk <= 1)
        print(f"DF toy a={a} slots={slots} profile={prof[:a+1]}: K_0=I {okI}, deg<=M(2L-1) {okdeg}, integrality on all discs {okInt}, "
              f"Hermite identity (D.9) entrywise exact {okH}, K_t rank-1 across groups {okloc}")
        loss = {t: -sum(prof[v2(s)] * (v2(s) - v2(t - s)) for s in S if s not in (0, t)) for t in c}
        print("    formula precision  1 + loss_t + (q_t-1)(a+1-i):", {t: 1 + loss[t] + (prof[v2(t)] - 1) * (a + 1 - v2(t)) for t in c},
              "  1+loss_t alone:", {t: 1 + loss[t] for t in c},
              "  min_t valuation of Hermite weights:", {t: min(v2(x) for x in c[t]) for t in c})
        # ---------------- (d) true minimal precision ----------------
        cols = [(t, j) for t in sorted(c) for j in range(prof[v2(t)])]
        rowsS = []; tgtS = []
        for (Y, X), P in K.items():
            jt = {t: jets(P, t, al, prof[v2(t)]) for t in c}
            rowsS.append([jt[t][j] for (t, j) in cols]); tgtS.append(peval(P, 0))
        # black-box class: deg <= N, integral on discs of radius 2^-(a+1) around every t in S: local saturation
        A = []
        for k in range(N + 1):
            mono = [Fr(0)] * k + [Fr(1)]
            A.append([x for t in S for x in (pcompose_affine(mono, t, al) + [Fr(0)] * (N + 1))[:N + 1]])
        # rows of A = images of monomials; lattice class = {c : cA integral}.  Basis: local row-reduce A to
        # unimodular * diag form and divide each pivot row by its pivot content.
        B = [r[:] for r in A]; m = len(B); n = len(B[0]); usedr = [False] * m; usedc = [False] * n
        Tm = [[Fr(int(i == j)) for j in range(m)] for i in range(m)]
        # saturate: repeatedly make the lattice spanned by rows 2-saturated (Gaussian elimination mod 2 on primitive rows)
        def prim(r):
            g = min(v2(x) for x in r if x != 0); return [x / Fr(2) ** g for x in r]
        B = [prim(r) for r in B]
        changed = True
        while changed:
            changed = False
            # find a nontrivial F_2-dependency among rows mod 2
            piv = {}; 
            for i in range(m):
                r = B[i][:]
                while True:
                    lead = next((c0 for c0 in range(n) if v2(r[c0]) == 0), None)
                    if lead is None: break
                    if lead in piv:
                        k = piv[lead]; f = r[lead] / B[k][lead]; r = [r[j] - f * B[k][j] for j in range(n)]
                    else: break
                if lead is None:
                    B[i] = prim(r); changed = True; break
                else:
                    B[i] = r; piv[lead] = i
        colidx = {(t, j): S.index(t) * (N + 1) + j for (t, j) in cols}
        rowsB = [[r[colidx[cj]] for cj in cols] for r in B]; tgtB = [r[S.index(0) * (N + 1)] for r in B]
        for label, rows_, tgt_ in (("black-box class (Lemma D.5 hypotheses)", rowsB, tgtB), ("structured class span{K_T(Y,X)}", rowsS, tgtS)):
            res = []
            for dlt in ([0] if a == 1 else [0, 4, 8, 12, 16, 24]):
                # precision p0 at odd nodes, p1 = p0 - dlt at level-1 nodes (scale level-1 columns by 2^dlt)
                cs = [dlt if v2(t) >= 1 else 0 for (t, j) in cols]
                need, rk = local_smith_minprec(rows_, tgt_, cs)
                res.append((dlt, need))
            print(f"    TRUE min precision, {label}: rank {rk}; (p0-p1, min p0) = {res}")
