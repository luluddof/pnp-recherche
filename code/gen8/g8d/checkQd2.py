"""Exhaustive check of the Q_{d+2} lemma: no set of 3 bad edges in Q_{d+2} is feasible; some 4-set is.
Two tau choices: random mod p, and tau_k = k+1 (special)."""
import itertools, sys
from lazy import Inst
for d in (1, 2, 3):
    l = d + 2
    for tag in ("random", "tau=k+1"):
        I = Inst(l, d)
        if tag != "random": I.tau = [k + 1 for k in range(l)]
        ne = len(I.E); allE = set(range(ne))
        f3 = sum(I.feasible(sorted(allE - set(B))) for B in itertools.combinations(range(ne), 3))
        f4 = 0
        for B in itertools.combinations(range(ne), 4):
            if I.feasible(sorted(allE - set(B))): f4 += 1
            if d == 3 and f4: break
        print(f"d={d} {tag}: feasible 3-bad sets {f3}, feasible 4-bad sets {f4}{'+' if d==3 else ''}", flush=True)
