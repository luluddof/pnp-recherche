#!/usr/bin/env python3
"""G18S r5: independent exact verification of G18e's explicit decompositions.
The term FUNCTIONS are taken from g18e/exact.py (two_block, copied verbatim into g18e_two_block.py); the VERIFIER is mine:
each term is tabulated on all 4^M end-bit assignments, flattened as (bits in group 0 | bits in group 1) of its direction,
and its exact rank over Q is computed (must be <= 1; with 2 groups per direction this is exactly 'factors over the groups');
the integer sum of the tables is compared entrywise with I."""
import itertools
from fractions import Fraction
from g18e_two_block import two_block
from r4_rank import rank_Q
def my_verify(name, u, wires, terms):
    M = len(wires); nb = 2 * M; party = [x for w in wires for x in w]
    allb = list(itertools.product((0, 1), repeat=nb)); tot = {b: 0 for b in allb}; maxrk = 0; dirs = []
    for (d, f) in terms:
        g0 = [k for k in range(nb) if (party[k] >> d) & 1 == 0]; g1 = [k for k in range(nb) if (party[k] >> d) & 1 == 1]
        Mat = []
        for a in itertools.product((0, 1), repeat=len(g0)):
            row = []
            for c in itertools.product((0, 1), repeat=len(g1)):
                b = [0] * nb
                for k, v in zip(g0, a): b[k] = v
                for k, v in zip(g1, c): b[k] = v
                v = f(tuple(b)); tot[tuple(b)] += v; row.append(v)
            Mat.append(row)
        maxrk = max(maxrk, rank_Q(Mat)); dirs.append(d)
    eq = all(tot[b] == int(all(b[2 * e] == b[2 * e + 1] for e in range(M))) for b in allb)
    print(f"{name}: M={M} N={len(terms)} directions used {sorted(set(dirs))} (counts {[dirs.count(x) for x in sorted(set(dirs))]})  max rank across groups = {maxrk}  sum == I over Z: {eq}")
    return eq and maxrk <= 1
if __name__ == "__main__":
    ok = True
    w = [(0, 3), (1, 2)];                 ok &= my_verify("u=2 slots{00,01}   ", 2, w, two_block(w, 2, 0, 1, [0], [1], []))
    w = [(0, 3), (3, 0), (1, 2)];         ok &= my_verify("u=2 slots{00,11,01}", 2, w, two_block(w, 2, 0, 1, [0], [2], [1]))
    w = [(0, 3), (3, 0), (1, 2), (2, 1)]; ok &= my_verify("u=2 full packet    ", 2, w, two_block(w, 2, 0, 1, [0, 1], [2, 3], []))
    w = [(0, 7), (1, 6), (2, 5)];         ok &= my_verify("u=3 three pairs    ", 3, w, two_block(w, 3, 2, 0, [0], [1], [2]))
    w = [(0, 7), (1, 6), (2, 5), (3, 4)]; ok &= my_verify("u=3 four pairs     ", 3, w, two_block(w, 3, 2, 0, [0, 2], [1, 3], []))
    print("ALL VERIFIED" if ok else "FAILURE")
