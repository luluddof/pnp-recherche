#!/usr/bin/env python3
"""Search-free version: the variable set is chosen from the HYPERGRAPH ONLY (signs unseen), as the
smallest set along (a) the natural order, (b) the min-degree peel order, whose induced sub-formula has
density >= RHO = 5.2 > ln2/ln(8/7) = 5.191.  With signs independent of the set, first moment gives
P(induced formula SAT) <= 2^v (7/8)^(RHO v) -> 0.  We report |S|/n and check UNSAT with a solver."""
import sys, json, math, time
from core_exp import gen, peel_order, solve
RHO = 5.2
def smallest_dense_prefix(cls, order, vmin=30):
    pos = {v: i for i, v in enumerate(order)}
    key = sorted(max(pos[abs(l)] for l in c) for c in cls)
    import bisect
    for p in range(vmin, len(order) + 1):
        if bisect.bisect_left(key, p) >= RHO * p: return p, bisect.bisect_left(key, p)
    return None, None
def induced(cls, order, p):
    S = set(order[:p]); return [c for c in cls if all(abs(l) in S for l in c)]
for n, D in [(400, 32), (400, 64), (400, 128), (1600, 64), (1600, 128), (1600, 256), (1600, 512), (6400, 256)]:
    out = []
    for sd in range(5):
        cls = gen(n, D, 7000 + sd)
        nat = list(range(1, n + 1)); po = peel_order(n, cls)
        p1, m1 = smallest_dense_prefix(cls, nat); p2, m2 = smallest_dense_prefix(cls, po)
        u1 = u2 = None
        if p2 <= 260: u2 = not solve(induced(cls, po, p2))[0]
        if p1 <= 260: u1 = not solve(induced(cls, nat, p1))[0]
        out.append((p1 / n, p2 / n, u1, u2))
    a = sorted(o[0] for o in out)[2]; b = sorted(o[1] for o in out)[2]
    print(json.dumps(dict(n=n, D=D, m=int(D*n), nat_frac=round(a, 4), sqrt_5p2_over_D=round(math.sqrt(RHO / D), 4),
          peel_frac=round(b, 4), peel_over_nat=round(b / a, 3), peel_times_D_pow=round(b * D ** 0.5, 3),
          unsat_checked_nat=[o[2] for o in out], unsat_checked_peel=[o[3] for o in out])), flush=True)
