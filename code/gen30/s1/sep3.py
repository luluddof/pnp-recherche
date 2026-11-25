import numpy as np, sys, time
from core import *
from sl import sl
k = 3; D = 64; V = np.diag(beta(k)); rng = np.random.default_rng(int(sys.argv[4]))
tests = [tuple(map(int, s.split(':'))) for s in sys.argv[2].split(',')]; starts = int(sys.argv[3]); mode = sys.argv[1]
mask = np.zeros((D, D))
if mode == 'bell': mask[0, :] = 1
else: mask[:] = 1
for q, r in tests:
    t0 = time.time(); best = 9e9; hit = None
    for s in range(starts):
        res, S, C = sl(V, q, r, mask, rng, 250, init='rand' if s else 'svd')
        best = min(best, res)
        if res < 1e-8:
            sup = np.argwhere(S); hit = (s, len(set(sup[:, 0])), int((sup[:, 0] == 0).sum())); break
    print(mode, f'q={q} r={r}: best res {best:.2e}', 'HIT start,#shifts,#Z' if hit else '', hit if hit else '', f'{time.time()-t0:.0f}s', flush=True)
