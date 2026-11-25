# k=4: from Majorana (10,46) add t free Paulis (Full or Bell); binary search min rank (heuristic). time-capped.
import numpy as np, sys, time
from core import *
from extend import run
from majorana import majorana_g
k = 4; D = 256; V = np.diag(beta(k)); g, supp, p0, J = majorana_g(k, 1, 2)
S0 = np.zeros((D, D), bool); S0[0, :] = supp; q0 = supp.sum(); c0 = np.zeros((D, D)); c0[0, :] = fwht(g) / D
mode = sys.argv[1]; ts = [int(x) for x in sys.argv[2].split(',')]; starts = int(sys.argv[3]); tcap = float(sys.argv[4])
fm = np.zeros((D, D))
if mode == 'bell': fm[0, :] = 1
else: fm[:] = 1
rng = np.random.default_rng(2); T0 = time.time()
print('base', q0, p0, flush=True)
for t in ts:
    lo, hi = 0, p0   # hi feasible
    while hi - lo > 1 and time.time() - T0 < tcap:
        mid = (lo + hi) // 2; ok = False
        for s in range(starts):
            res, keep = run(V, S0, t, mid, fm, rng, its=120, init=c0 if s == 0 else None)
            if res < 1e-8: ok = True; break
        if ok: hi = mid
        else: lo = mid
        print(f'   {mode} t={t} r={mid} ok={ok} res={res:.2e} ({time.time()-T0:.0f}s)', flush=True)
    print(f'{mode}: q={q0}+{t}={q0+t}: min rank found {hi}', flush=True)
