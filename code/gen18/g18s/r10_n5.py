#!/usr/bin/env python3
"""G18S r10: the 5-term mod-2 DF identity for 'u=3 four pairs' (wires (0,7),(1,6),(2,5),(3,4)), directions (0,0,1,1,2).
Re-find it, verify it mod 2 with an independent tabulation, save it, and test whether it lifts:
 (a) is sum_k tau_k == I over Z ?  (b) is there kappa in {+-1}^5 (or kappa in Q, by exact linear solve) with sum kappa_k tau_k == I over Z ?"""
import itertools, json
from fractions import Fraction
from dflib import *
u = 3; wires = [(0, 7), (1, 6), (2, 5), (3, 4)]; dirs = (0, 0, 1, 1, 2)
r, tabs = sat_search(u, wires, dirs); assert r == 'SAT'
M = len(wires); nb = 2 * M; allb = list(itertools.product((0, 1), repeat=nb))
I = [int(all(b[2 * e] == b[2 * e + 1] for e in range(M))) for b in allb]
T = []
for k, d in enumerate(dirs):
    g0, g1 = groups(u, wires, d)
    T.append([tabs[('t', k, 0, tuple(b[i] for i in g0))] * tabs[('t', k, 1, tuple(b[i] for i in g1))] for b in allb])
print("mod 2:", all(sum(T[k][j] for k in range(5)) % 2 == I[j] for j in range(len(allb))))
S = [sum(T[k][j] for k in range(5)) for j in range(len(allb))]
print("over Z with kappa = 1:", S == I, " max entry of sum:", max(S), " #entries where sum != I:", sum(a != b for a, b in zip(S, I)))
print("supports (number of nonzero entries) of the 5 terms:", [sum(t) for t in T])
best = None
for sg in itertools.product((1, -1), repeat=5):
    if [sum(s * T[k][j] for k, s in enumerate(sg)) for j in range(len(allb))] == I: best = sg; break
print("signs kappa in {+-1}^5 giving an exact identity over Z:", best)
json.dump({"u": u, "wires": wires, "dirs": dirs, "tables": {str(k): v for k, v in tabs.items()}}, open("n5_gf2_solution.json", "w"))
# number of solutions with integer sum == I: resample a few SAT models? (single model only)
