#!/usr/bin/env python3
"""G17a part 4: TRUE minimal precision for recovering F(0) mod 2 from Taylor data, by exact local
(2-adic) linear algebra, for (BB) the black-box class of Lemma D.5 and (ST) the Z_(2)-span of the
kernel entries K_T(Y,X) of one packet.  Compared with: formula (D.9), the prover's 'loss', and the
product-witness bound g_t."""
from fractions import Fraction as Fr
import sys, itertools
from rev1_lib import *
from rev2 import kernel, jets, local_smith_minprec

def saturated_basis(A):
    R = [r[:] for r in A]; m = len(R); n = len(R[0]); used = [False] * m; basis = []
    for _ in range(m):
        best = None
        for i in range(m):
            if used[i]: continue
            for c in range(n):
                if R[i][c] != 0:
                    v = v2(R[i][c])
                    if best is None or v < best[0]: best = (v, i, c)
        if best is None: break
        v, i, c = best; used[i] = True
        for k in range(m):
            if used[k] or R[k][c] == 0: continue
            f = R[k][c] / R[i][c]
            R[k] = [R[k][j] - f * R[i][j] for j in range(n)]
        basis.append([x / Fr(2) ** v for x in R[i]])
    return basis

def run(a, slots, prof, deltas):
    L, S = grid(a); al = 2 ** (a + 1); M = len(slots); N = M * (2 * L - 1)
    nodes = [t for t in S if t != 0 and prof[v2(t)] > 0]
    cols = [(t, j) for t in nodes for j in range(prof[v2(t)])]
    K = kernel(a, slots)
    rowsS = []; tgtS = []
    for (Y, X), P in K.items():
        jt = {t: jets(P, t, al, prof[v2(t)]) for t in nodes}
        rowsS.append([jt[t][j] for (t, j) in cols]); tgtS.append(peval(P, 0))
    Wd = max(N + 1, max(prof)); A = []
    for k in range(N + 1):
        mono = [Fr(0)] * k + [Fr(1)]
        A.append([x for t in S for x in (pcompose_affine(mono, t, al) + [Fr(0)] * Wd)[:Wd]])
    B = saturated_basis(A)
    ci = {(t, j): S.index(t) * Wd + j for (t, j) in cols}
    rowsB = [[r[ci[c]] for c in cols] for r in B]; tgtB = [r[S.index(0) * Wd] for r in B]
    loss = {i: -sum(prof[v2(s)] * (v2(s) - v2(t - s)) for s in S if s not in (0, t)) for t in nodes for i in [v2(t)]}
    print(f"a={a} L={L} wires={M} N={N} profile={prof[:a+1]} (sum q = {len(cols)}):")
    print("   formula (D.9) needs p_i >= 1+loss_i+(q_i-1)(a+1-i):", {i: 1 + loss[i] + (prof[i] - 1) * (a + 1 - i) for i in sorted(loss)},
          "| prover's VT hypothesis p_i >= loss_i:", {i: loss[i] for i in sorted(loss)})
    for label, rows_, tgt_ in (("BB", rowsB, tgtB), ("ST", rowsS, tgtS)):
        res = []
        for d in deltas:
            cs = [d if v2(t) >= 1 else 0 for (t, j) in cols]
            need, rk = local_smith_minprec(rows_, tgt_, cs)
            res.append((need, need - d))
        print(f"   {label}: exact minimal (p_0, p_1) pairs along p_0-p_1 = {deltas}: {res}")
    sys.stdout.flush()

if __name__ == "__main__":
    run(1, [0, 1], [6, 0, 0], [0])
    run(1, [0, 1], [4, 0, 0], [0])
    run(1, [0, 1], [40, 0, 0], [0])          # heavy over-determination: formula loss grows, true need must not
    run(2, [0, 1], [3, 2, 0, 0], [0, 2, 4, 6, 8, 12])
    run(2, [0, 1], [6, 6, 0, 0], [0, 2, 4, 6, 8, 12])
    run(2, [0, 1, 2], [9, 5, 0, 0], [0, 4, 8, 12, 16, 20])
    run(2, [0, 1, 3], [9, 5, 0, 0], [0, 4, 8, 12, 16, 20])
