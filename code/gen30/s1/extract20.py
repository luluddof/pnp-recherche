import numpy as np, pickle
from core import *
from extend import run
from majorana import majorana_g
k = 3; D = 64; V = np.diag(beta(k)); g, supp, p0, J = majorana_g(k, 1, 1)
S0 = np.zeros((D, D), bool); S0[0, :] = supp
fm = np.ones((D, D)); rng = np.random.default_rng(11); c0 = np.zeros((D, D)); c0[0, :] = fwht(g) / D
for r in (6, 5, 4):
  for t in (12, 16, 20, 24):
    found = False
    for s in range(6):
        res, keep = run(V, S0, t, r, fm, rng, init=c0 if s == 0 else None)
        if res < 1e-8: found = True; break
    print('r', r, 't', t, 'found', found, 'res', res, flush=True)
    if found:
        # recompute C on support keep with final L
        pickle.dump(keep, open(f'keep_r{r}_t{t}.pkl', 'wb'))
        sup = np.argwhere(keep); shifts = sorted(set(sup[:, 0])); print('  shifts', shifts, 'count per shift', [int((sup[:, 0] == a).sum()) for a in shifts])
        break
