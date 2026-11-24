#!/usr/bin/env python3
"""G18S r6: exact optimisation of the flattening bound on the uniform 2-symbol instance U(u,2)
(labels {0,1}^u, one wire per ordered (row label, column label) pair, M = 4^u, direction i groups by bit i).
With 2 groups per direction:  w - c_i = |sum over selected wires of (bit_i(row-index end) - bit_i(col-index end))|.
 (a) vertex cuts by label classes (the notes' brute force);
 (b) vertex cuts where the row vertices and the column vertices of a label class are assigned separately;
 (c) ALL oriented wire selections (class level, integer): exact DP over reachable sum vectors;
 (d) LP relaxation of (c) (fractional selections), u <= 5, via scipy."""
import itertools, numpy as np, math
from scipy.optimize import linprog
def labels(u): return list(itertools.product((0, 1), repeat=u))
for u in (1, 2, 3):
    L = labels(u); n = len(L); M = 4 ** u
    # (a)
    best = 0
    for asg in itertools.product((0, 1, 2), repeat=n):
        A = [L[k] for k in range(n) if asg[k] == 1]; B = [L[k] for k in range(n) if asg[k] == 2]
        # wires A->B oriented (row idx end in A): both (row in A, col in B) and (row in B, col in A): 2|A||B| wires
        val = min(abs(2 * (len(B) * sum(x[i] for x in A) - len(A) * sum(y[i] for y in B))) for i in range(u)) if A and B else 0
        best = max(best, val)
    print(f"u={u}: (a) label-class vertex cuts: max = {best}/{M} = {best/M:.5f}")
    # (b) stats per side assignment: (n1, S1 vector, n2, S2 vector)
    X = np.array(L)
    asgs = np.array(list(itertools.product((0, 1, 2), repeat=n)), dtype=np.int64)
    n1 = (asgs == 1).sum(1); n2 = (asgs == 2).sum(1); S1 = (asgs == 1) @ X; S2 = (asgs == 2) @ X
    bestb = 0
    for k in range(len(asgs)):   # row side assignment k, all column side assignments vectorised
        # delta=+1: row in R1, col in C2 ; delta=-1: row in R2, col in C1
        f = S1[k][None, :] * n2[:, None] - n1[k] * S2 - (S2[k][None, :] * n1[:, None] - n2[k] * S1)
        bestb = max(bestb, int(np.abs(f).min(1).max()))
    print(f"u={u}: (b) row/column-separate label-class vertex cuts: max = {bestb}/{M} = {bestb/M:.5f}")
    # (c) DP over reachable sum vectors, classes (x,y) with vector v = x - y
    reach = {tuple([0] * u)}
    for x in L:
        for y in L:
            v = tuple(a - b for a, b in zip(x, y))
            if not any(v): continue
            new = set(reach)
            for s in reach:
                new.add(tuple(a + b for a, b in zip(s, v))); new.add(tuple(a - b for a, b in zip(s, v)))
            reach = new
    bestc = max(min(abs(a) for a in s) for s in reach)
    print(f"u={u}: (c) all oriented class selections: max = {bestc}/{M} = {bestc/M:.5f}   [C(2u,u)/4^u = {math.comb(2*u,u)/4**u:.5f}]")
for u in (1, 2, 3, 4, 5):
    L = labels(u); V = np.array([[a - b for a, b in zip(x, y)] for x in L for y in L], dtype=float)  # (4^u, u)
    nc = len(V)
    # maximise t s.t. t <= sum_c delta_c V[c,i] for all i, -1 <= delta <= 1  (all-positive sign pattern; others equivalent by bit flips)
    c = np.zeros(nc + 1); c[-1] = -1
    A = np.hstack([-V.T, np.ones((u, 1))]); b = np.zeros(u)
    res = linprog(c, A_ub=A, b_ub=b, bounds=[(-1, 1)] * nc + [(None, None)], method="highs")
    print(f"u={u}: (d) LP optimum = {-res.fun:.6f}/{4**u} = {-res.fun/4**u:.6f}   C(2u,u)/4^u = {math.comb(2*u,u)/4**u:.6f}")
