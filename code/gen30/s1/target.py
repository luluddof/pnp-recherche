# test feasibility of rank r with q Paulis (Full or Bell), several starts; prints best residual
import numpy as np, sys, time
from core import *
from sl import sl
k = int(sys.argv[1]); mode = sys.argv[2]; pairs = [tuple(map(int, s.split(':'))) for s in sys.argv[3].split(',')]; starts = int(sys.argv[4])
its = int(sys.argv[5]) if len(sys.argv) > 5 else 300
D = 4**k; V = np.diag(beta(k)); rng = np.random.default_rng(int(time.time()) % 1000)
mask = np.zeros((D, D)); 
if mode == 'bell': mask[0, :] = 1
else: mask[:] = 1
for q, r in pairs:
    t0 = time.time(); best = 1e9; nfound = 0
    for s in range(starts):
        res, S, C = sl(V, q, r, mask, rng, its, init='rand' if s else 'svd')
        best = min(best, res)
        if res < 1e-8:
            nfound += 1; sup = np.argwhere(S); nsh = len(set(sup[:, 0]))
            print(f'  FOUND q={q} r={r} start {s}: #X-shifts {nsh}, shifts {sorted(set(sup[:,0]))[:12]}', flush=True); break
    print(f'{mode} k={k} q={q} r={r}: best residual {best:.3e} ({time.time()-t0:.0f}s)', flush=True)
