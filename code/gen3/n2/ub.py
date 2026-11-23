"""Upper bound on the FULL cover-LP (all found bad covers): v_C = 1/(min size of a bad cover containing C) is LP-feasible."""
import sys, random, numpy as np, fko
from collections import defaultdict
for n, a in [(50,10),(100,7.4),(100,12),(200,12),(300,7.4)]:
    m = int(a*n**1.5); rng = random.Random(1000*n+int(10*a)); V, S = fko.gen(n, m, rng)
    B12, B3 = fko.spectral(n, V, S) if n <= 200 else (None, None)
    covs, _ = fko.cycle_covers(n, V, 30 if n <= 100 else 8, rng)
    sgn = np.prod(S, axis=1); ms = {}
    for K in covs:
        if np.prod([sgn[k] for k in K]) == -1:
            for k in K: ms[k] = min(ms.get(k, 99), len(K))
    ub = sum(1/v for v in ms.values())
    print(f"n={n} a15={a} m={m} LPfull<= {ub:.0f}  need=B12/2={B12/2 if B12 else float('nan'):.0f}  spec_ok={B12+B3<m if B12 else None}", flush=True)
