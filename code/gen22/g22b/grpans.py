# Group-invariant ansatz: rows/cols = group Gm (order N), cell (a, a*c) gets colour (B(c), a*f(c)); invariant under left mult.
import sys, itertools, time, json
from pysat.solvers import Cadical153
from pysat.card import CardEnc
from shape import verify_shape, show

def groups(N):
    import sympy.combinatorics.named_groups as ng
    from sympy.combinatorics import PermutationGroup
    out = {}
    if N == 8:
        from sympy.combinatorics.group_constructs import DirectProduct
        out['Z8'] = ng.CyclicGroup(8); out['Z4xZ2'] = DirectProduct(ng.CyclicGroup(4), ng.CyclicGroup(2)); out['D4'] = ng.DihedralGroup(4)
        from sympy.combinatorics import Permutation
        # quaternion group Q8 as permutations of 8 points (regular rep)
        i = Permutation([2,3,1,0,6,7,5,4]); j = Permutation([4,5,7,6,1,0,2,3])
        Q = PermutationGroup([i, j])
        if Q.order() == 8: out['Q8'] = Q
    return out
def table(Gp):
    el = list(Gp.generate()); idx = {e: k for k, e in enumerate(el)}
    e0 = idx[Gp.identity]; el = [el[e0]] + [x for k, x in enumerate(el) if k != e0]; idx = {e: k for k, e in enumerate(el)}
    mul = [[idx[el[a] * el[b]] for b in range(len(el))] for a in range(len(el))]
    inv = [idx[x**-1] for x in el]
    return mul, inv
def solve(mul, inv, nG, nH, size):
    N = len(mul); nb = nG + nH; vid = itertools.count(1)
    F = [[next(vid) for v in range(N)] for c in range(N)]; B = [[next(vid) for i in range(nb)] for c in range(N)]
    cl = []
    for c in range(N):
        cl.append(F[c][:]); cl.append(B[c][:])
        for v, w in itertools.combinations(range(N), 2): cl.append([-F[c][v], -F[c][w]])
        for i, j in itertools.combinations(range(nb), 2): cl.append([-B[c][i], -B[c][j]])
    top = next(vid)
    for i in range(nb):
        enc = CardEnc.equals(lits=[B[c][i] for c in range(N)], bound=size, top_id=top); top = max(top, enc.nv) + 1; cl += enc.clauses
    vid = itertools.count(top + 1); SB = {}
    for c1, c2 in itertools.permutations(range(N), 2):
        if (c2, c1) in SB: SB[(c1, c2)] = SB[(c2, c1)]; continue
        s = next(vid); SB[(c1, c2)] = s
        for i in range(nb): cl.append([-B[c1][i], -B[c2][i], s]); cl.append([-s, -B[c1][i], B[c2][i]])
    for c1, c2 in itertools.permutations(range(N), 2):
        s = SB[(c1, c2)]
        for v1 in range(N):
            for v2 in range(N):
                pre = [-s, -F[c1][v1], -F[c2][v2]]
                if v1 == v2: cl.append(pre); continue
                if mul[inv[v1]][c1] == mul[inv[v2]][c2]: cl.append(pre); continue
                c = mul[mul[v1][inv[v2]]][c2]; cp = mul[mul[v2][inv[v1]]][c1]
                cl.append(pre + [SB[(c, cp)]] if c != cp else pre)
                for i in range(nb):
                    for j in range(nb):
                        if (i < nG) == (j < nG): cl.append(pre + [-B[c1][i], -B[c][j]])
                for u in range(N):   # v1^-1 f(c) = v2^-1 f(cp)  <=> f(cp) = v2 v1^-1 f(c)
                    cl.append(pre + [-F[c][u], F[cp][mul[mul[v2][inv[v1]]][u]]])
    S = Cadical153(bootstrap_with=cl)
    return S, F, B
if __name__ == "__main__":
    N, nG, nH, size = map(int, sys.argv[1:5])
    for name, Gp in groups(N).items():
        mul, inv = table(Gp); S, F, B = solve(mul, inv, nG, nH, size)
        found = None
        for t in range(200):
            if not S.solve(): break
            m = set(l for l in S.get_model() if l > 0)
            f = [next(v for v in range(N) if F[c][v] in m) for c in range(N)]; bl = [next(i for i in range(nG+nH) if B[c][i] in m) for c in range(N)]
            col = [[None]*N for _ in range(N)]
            for a in range(N):
                for c in range(N):
                    b = mul[a][c]; col[a][b] = ('S' if bl[c] < nG else 'T', bl[c]*N + mul[a][f[c]])
            from transl import sign_solve
            # reuse the F_2 sign solver on an explicit colouring
            import transl
            cells = transl_sign(col, N) if False else None
            # inline sign solve
            rows = []
            for a, a2 in itertools.combinations(range(N), 2):
                for b in range(N):
                    for b2 in range(N):
                        if b != b2 and col[a][b] == col[a2][b2]:
                            v = 0
                            for (x, y) in ((a, b), (a2, b2), (a, b2), (a2, b)): v ^= 1 << (x*N + y)
                            rows.append((v, 1))
            piv = {}; bad = False
            for v, r in rows:
                while v:
                    h = v.bit_length() - 1
                    if h in piv: pv, pr = piv[h]; v ^= pv; r ^= pr
                    else: piv[h] = (v, r); break
                if v == 0 and r == 1: bad = True; break
            if not bad:
                x = 0
                for h in sorted(piv):
                    pv, pr = piv[h]; val = pr ^ (bin(pv & x & ~(1 << h)).count('1') & 1)
                    if val: x |= 1 << h
                found = [[(col[a][b][0], col[a][b][1], -1 if (x >> (a*N+b)) & 1 else 1) for b in range(N)] for a in range(N)]; break
            S.add_clause([-F[c][f[c]] for c in range(N)] + [-B[c][bl[c]] for c in range(N)])
        print(name, "colourings tried", t, "FOUND" if found else "none", flush=True)
        if found:
            show(N, found); verify_shape(N, found); json.dump({"N": N, "group": name, "cells": found}, open(f"grp_{N}_{name}.json", "w"))
