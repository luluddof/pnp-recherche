import numpy as np, pickle, time
from core import *
from sl import trunc
from extend import run
from majorana import majorana_g
k = 4; D = 256; V = np.diag(beta(k)); g, supp, p0, J = majorana_g(k, 1, 2)
S0 = np.zeros((D, D), bool); S0[0, :] = supp; c0 = np.zeros((D, D)); c0[0, :] = fwht(g) / D
rng = np.random.default_rng(2); T0 = time.time()
full = np.ones((D, D)); found = None
for s in range(6):
    res, keep = run(V, S0, 24, 38, full, rng, its=120, init=c0 if s == 0 else None)
    print('full r=38 start', s, 'res', res, f'{time.time()-T0:.0f}s', flush=True)
    if res < 1e-8: found = keep; break
if found is not None:
    pickle.dump(found, open('k4_full_q34_r38.pkl', 'wb'))
    sup = np.argwhere(found); nb = int((sup[:, 0] != 0).sum())
    print('support', len(sup), 'non-Bell terms', nb, 'shifts', sorted(set(int(a) for a in sup[:, 0])))
    # re-polish from random inits on the fixed support, several ranks
    for r in (38, 37, 36):
        ok = False
        for s in range(3):
            L = trunc(rng.standard_normal((D, D)), r)
            for it in range(6000):
                C = pinv_trans(ptrans(V - L) * found); L = trunc(V - C, r)
                if it % 500 == 499:
                    res = np.linalg.norm(V - C - L)
                    if res < 1e-10: break
            res = np.linalg.norm(V - C - L)
            if res < 1e-8:
                ok = True; c = ptrans(C); offn = np.abs(c[1:]).max(); print(f'  fixed support r={r}: feasible; max |non-Bell coef| = {offn:.3e}', flush=True); break
        if not ok: print(f'  fixed support r={r}: not reached (res {res:.2e})', flush=True)
