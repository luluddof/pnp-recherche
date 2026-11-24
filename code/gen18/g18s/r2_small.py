#!/usr/bin/env python3
"""G18S r2: the notes' small cases.  Bounds by brute force over ALL vertex cuts (A') and ALL oriented wire subsets (A''),
then exact GF(2) search (SAT, p = 1, kappa = 1, 0/1 tables) for N below / at the bound."""
import itertools, sys, time, multiprocessing as mp
from dflib import *
INST = [("u=1 both slots", 1, [(0, 1), (1, 0)]),
        ("u=2 slots{00,01}", 2, [(0, 3), (1, 2)]),
        ("u=2 slots{00,11}", 2, [(0, 3), (3, 0)]),
        ("u=2 slots{00,11,01}", 2, [(0, 3), (3, 0), (1, 2)]),
        ("u=2 full packet", 2, [(0, 3), (3, 0), (1, 2), (2, 1)]),
        ("u=3 three pairs", 3, [(0, 7), (1, 6), (2, 5)]),
        ("u=3 four pairs", 3, [(0, 7), (1, 6), (2, 5), (3, 4)])]
def job(args):
    u, wires, dirs = args
    t = time.time(); r, tabs = sat_search(u, wires, dirs)
    ok = verify_gf2(u, wires, dirs, tabs) if r == 'SAT' else None
    return dirs, r, ok, time.time() - t
if __name__ == "__main__":
    for name, u, wires in INST:
        b = bounds(u, wires)
        print(f"{name}: M={len(wires)}  log2 N >=  ThmA {b['thmA']}  A'(vertex cuts) {b['Ap']}  A''(oriented subsets) {b['App']}   arg A'' {b['argApp']}", flush=True)
