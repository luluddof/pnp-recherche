# From a Bell support S0 (Majorana), add t free Paulis (Full: any; Bell: Z-type) and try to reach rank r.
import numpy as np, sys, time
from core import *
from sl import trunc
from majorana import majorana_g
def run(V, S0, t, r, freemask, rng, its=200, init=None):
    D = V.shape[0]
    L = trunc(V - pinv_trans(init), r) if init is not None and r > 0 else trunc(rng.standard_normal((D, D)), r)
    for it in range(its):
        c = ptrans(V - L); keep = S0.copy()
        if t > 0:
            cf = np.abs(c) * freemask * (~S0); idx = np.argsort(-cf.ravel())[:t]
            keep.ravel()[idx] = True
        C = pinv_trans(c * keep); L = trunc(V - C, r)
    hist = []
    for it in range(4000):
        C = pinv_trans(ptrans(V - L) * keep); L = trunc(V - C, r)
        res = np.linalg.norm(V - C - L); hist.append(res)
        if res < 1e-11: break
        if it > 300 and res > 0.995 * hist[-200] and res > 1e-4: break
    return res, keep
if __name__ == '__main__':
    k, d, w0 = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); ts = [int(x) for x in sys.argv[4].split(',')]; starts = int(sys.argv[5])
    D = 4**k; V = np.diag(beta(k)); g, supp, p0, J = majorana_g(k, d, w0)
    S0 = np.zeros((D, D), bool); S0[0, :] = supp; q0 = supp.sum()
    c0 = np.zeros((D, D)); c0[0, :] = fwht(g) / D
    print('base Majorana q0', q0, 'p0', p0, flush=True)
    rng = np.random.default_rng(7)
    for mode in ('bell', 'full'):
        fm = np.zeros((D, D)); 
        if mode == 'bell': fm[0, :] = 1
        else: fm[:] = 1
        for t in ts:
            best = p0
            for r in range(p0 - 1, -1, -1):
                ok = False
                for s in range(starts):
                    res, keep = run(V, S0, t, r, fm, rng, init=c0 if s == 0 else None)
                    if res < 1e-8: ok = True; break
                if not ok: break
                best = r; nsh = len(set(np.argwhere(keep)[:, 0]))
            print(f'{mode}: q = {q0}+{t} = {q0+t}: min rank found {best}' + (f' (#X-shifts {nsh})' if best < p0 else ''), flush=True)
