"""Start from the layered good set (piecewise construction) and greedily add bad edges keeping
generic consistency (mod p, random tau). Reports whether layered sets are maximal."""
import sys, random, numpy as np
from core import edges
from greedy_max import Ech, p
from math import comb
def run(l, d, r, seed=0, trials=3):
    E = edges(l); n = (1 << l) * d; rng = random.Random(seed)
    res = []
    for tr in range(trials):
        taus = [rng.randrange(2, p - 1) for _ in range(l)]
        def row(z, w, k):
            rr = np.zeros(n + 1, dtype=np.int64); pw = taus[k]
            for j in range(d):
                rr[z*d+j] = (rr[z*d+j] + pw) % p; rr[w*d+j] = (rr[w*d+j] + pw) % p; pw = pw * taus[k] % p
            rr[n] = p - 2; return rr
        bl = lambda w: (w - r - 1) // (d + 1)
        good = [e for e in E if bl(bin(e[0]).count('1')) == bl(bin(e[0]).count('1') + 1)]
        bad = [e for e in E if bl(bin(e[0]).count('1')) != bl(bin(e[0]).count('1') + 1)]
        ech = Ech(n)
        for e in good:
            ok, _ = ech.add(row(*e)); assert ok
        rng.shuffle(bad); added = 0
        for e in bad:
            ok, _ = ech.add(row(*e)); added += ok
        res.append(added)
    return len(bad), res
for l in range(5, int(sys.argv[1]) + 1):
    for d in (2, 3):
        best = min(range(d + 1), key=lambda r: sum(comb(l - 1, w) for w in range(l) if (w - r - 1) // (d + 1) != (w - r) // (d + 1)))
        nb, res = run(l, d, best)
        print(f"l={l} d={d} r={best}: layered bad total={nb}; extra edges addable={res}", flush=True)
