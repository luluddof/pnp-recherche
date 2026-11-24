#!/usr/bin/env python3
"""G17a part 5: structured vs black-box class with several packets (a=1, two full packets, M=4)
and one full packet at a=2 (M=4)."""
from fractions import Fraction as Fr
import sys
from rev1_lib import *
from rev2 import kernel, jets, local_smith_minprec
from rev4 import saturated_basis
import rev4
def run_multi(a, packets, prof, deltas):
    L, S = grid(a); al = 2 ** (a + 1); M = sum(len(p) for p in packets); N = M * (2 * L - 1)
    polys = [[Fr(1)]]
    for p in packets:
        K = kernel(a, p); new = []
        seen = set()
        for P in K.values():
            key = tuple(P)
            if key in seen: continue
            seen.add(key)
            for Q in polys: new.append(pmul(P, Q))
        # dedupe
        d = {}
        for P in new: d[tuple(P)] = P
        polys = list(d.values())
    nodes = [t for t in S if t != 0 and prof[v2(t)] > 0]
    cols = [(t, j) for t in nodes for j in range(prof[v2(t)])]
    rowsS = []; tgtS = []
    for P in polys:
        jt = {t: jets(P, t, al, prof[v2(t)]) for t in nodes}
        rowsS.append([jt[t][j] for (t, j) in cols]); tgtS.append(peval(P, 0))
    Wd = max(N + 1, max(prof)); A = []
    for k in range(N + 1):
        mono = [Fr(0)] * k + [Fr(1)]
        A.append([x for t in S for x in (pcompose_affine(mono, t, al) + [Fr(0)] * Wd)[:Wd]])
    B = saturated_basis(A)
    ci = {(t, j): S.index(t) * Wd + j for (t, j) in cols}
    rowsB = [[r[ci[c]] for c in cols] for r in B]; tgtB = [r[S.index(0) * Wd] for r in B]
    print(f"a={a} packets={packets} M={M} N={N} profile={prof[:a+1]}  ({len(polys)} distinct kernel polynomials)")
    for label, rows_, tgt_ in (("BB", rowsB, tgtB), ("ST", rowsS, tgtS)):
        res = []
        for d in deltas:
            cs = [d if v2(t) >= 1 else 0 for (t, j) in cols]
            need, rk = local_smith_minprec(rows_, tgt_, cs); res.append((need, need - d))
        print(f"   {label}: rank {rk}; exact minimal (p_0,p_1) along p_0-p_1={deltas}: {res}"); sys.stdout.flush()
run_multi(1, [[0, 1], [0, 1]], [8, 0, 0], [0])
run_multi(1, [[0, 1], [0, 1], [0, 1]], [10, 0, 0], [0])
run_multi(2, [[0, 1, 2, 3]], [8, 4, 0, 0], [0, 6, 12, 18])
