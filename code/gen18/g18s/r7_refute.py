#!/usr/bin/env python3
"""G18S r7: refutation attempt for the generalised bound A'' where it is STRICTLY stronger than every vertex cut.
Enumerate small instances (u=2: M<=5; u=3: M<=4, wires between arbitrary corners), keep those with A'' > A',
and run the exact GF(2) search for N = 2^{A''} - 1 terms over all direction profiles (a SAT answer would refute A'')."""
import itertools, random, sys, time
from dflib import *
random.seed(7); seen = set(); cases = []
for u, Ms in ((2, (2, 3, 4, 5)), (3, (3, 4))):
    P = 2 ** u; pairs = [(p, q) for p in range(P) for q in range(P) if p != q]
    for M in Ms:
        for _ in range(300):
            wires = tuple(sorted(random.choice(pairs) for _ in range(M)))
            if (u, wires) in seen: continue
            seen.add((u, wires)); b = bounds(u, list(wires))
            if b['App'] > b['Ap']: cases.append((u, wires, b))
print(f"{len(seen)} instances, {len(cases)} with A'' > A' (vertex cuts)")
cases.sort(key=lambda c: (-c[2]['App'], len(c[1])))
t0 = time.time(); done = 0; sat = 0
for u, wires, b in cases:
    if time.time() - t0 > 400: break
    N = 2 ** b['App'] - 1
    if N > 3: continue
    rs = []
    for dirs in itertools.combinations_with_replacement(range(u), N):
        r, tabs = sat_search(u, list(wires), dirs); rs.append(r)
        if r == 'SAT': sat += 1; print("REFUTATION CANDIDATE", u, wires, dirs, verify_gf2(u, list(wires), dirs, tabs))
    done += 1
    if done <= 12: print(f"u={u} wires={wires}  ThmA {b['thmA']} A' {b['Ap']} A'' {b['App']}  N={N}: {sorted(set(rs))}", flush=True)
print(f"tested {done} instances at N = 2^A'' - 1: SAT found in {sat}")
