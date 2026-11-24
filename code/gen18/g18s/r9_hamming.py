#!/usr/bin/env python3
"""G18S r9: loads concentrated on Hamming-distance-1 pair classes: DF value is exactly 2^{min_i M_i} (<= 2^{M/u})."""
import itertools
from dflib import *
for name, u, wires in [("u=2, 2 wires per direction", 2, [(0, 1), (2, 3), (0, 2), (1, 3)]),
                       ("u=3, 1 wire per direction", 3, [(0, 1), (0, 2), (0, 4)]),
                       ("u=3, 1 wire per direction (spread)", 3, [(0, 1), (2, 0), (7, 3)]),
                       ("u=2, M=(3,2)", 2, [(0, 1), (2, 3), (3, 2), (0, 2), (1, 3)])]:
    b = bounds(u, wires); Mi = [sum(1 for (p, q) in wires if (p ^ q) >> d & 1) for d in range(u)]
    lb = 2 ** b['App']; out = []
    for N in (lb - 1, lb):
        if N < 1: continue
        r = 'UNSAT'
        for dirs in itertools.combinations_with_replacement(range(u), N):
            rr, t = sat_search(u, wires, dirs)
            if rr == 'SAT': r = f"SAT{dirs} verified={verify_gf2(u, wires, dirs, t)}"; break
        out.append(f"N={N}: {r}")
    print(f"{name}: M={len(wires)} crossing counts M_i={Mi}  ThmA {b['thmA']} A' {b['Ap']} A'' {b['App']} (= min M_i: {b['App']==min(Mi)})  GF(2): " + "; ".join(out))
