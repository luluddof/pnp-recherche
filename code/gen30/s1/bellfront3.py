import numpy as np, time, sys
from core import *
from sl import sl
k = 3; D = 64; V = np.diag(beta(k)); rng = np.random.default_rng(3)
mode = sys.argv[1]
mask = np.zeros((D, D))
if mode == 'bell': mask[0, :] = 1
else: mask[:] = 1
qs = [int(x) for x in sys.argv[2].split(',')]; prev = int(sys.argv[3]); starts = int(sys.argv[4])
for q in qs:
    best = None; r = prev
    while r >= 0:
        ok = False
        for s in range(starts):
            res, S, C = sl(V, q, r, mask, rng, 150, init='rand' if s else 'svd')
            if res < 1e-8: ok = True; break
        if not ok: break
        best = r; nsh = len(set(np.argwhere(S)[:, 0])); r -= 1
    if best is not None: prev = best
    print(mode, 'k=3 q', q, 'min r found', best, 'shifts in last', nsh if best is not None else '-', flush=True)
