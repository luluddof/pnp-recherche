# sparse (Pauli) + low-rank alternating projections: V = C + L, |supp C| <= q (within mask), rank L <= r
import numpy as np
from core import *
def trunc(M, r):
    if r == 0: return np.zeros_like(M)
    U, s, Vt = np.linalg.svd(M); return (U[:, :r] * s[:r]) @ Vt[:r]
def sl(V, q, r, mask, rng, its=400, init='rand'):
    D = V.shape[0]
    if init == 'rand': L = trunc(rng.standard_normal((D, D)), r) * np.linalg.norm(V) / np.sqrt(D)
    else: L = trunc(V, r)
    for t in range(its):
        c = ptrans(V - L) * mask
        idx = np.argsort(-np.abs(c).ravel())[:q]; ck = np.zeros_like(c).ravel(); ck[idx] = c.ravel()[idx]; ck = ck.reshape(c.shape)
        C = pinv_trans(ck); L = trunc(V - C, r)
    # polish with fixed support
    S = np.abs(ck) > 0
    hist = []
    for t in range(4000):
        C = pinv_trans(ptrans(V - L) * S); L = trunc(V - C, r)
        res = np.linalg.norm(V - C - L); hist.append(res)
        if res < 1e-11: break
        if t > 300 and res > 0.995 * hist[-200] and res > 1e-4: break
    return res, S, C
def frontier(k, bell, rmax, qlist, starts, seed=0, its=300):
    rng = np.random.default_rng(seed); D = 4**k; V = np.diag(beta(k))
    mask = np.zeros((D, D)); 
    if bell: mask[0, :] = 1
    else: mask[:] = 1
    out = {}
    for q in qlist:
        best = None
        for r in range(max(0, int(np.ceil(D/(4*q)))), rmax+1):
            ok = False
            for s in range(starts):
                res, S, C = sl(V, q, r, mask, rng, its, init='rand' if s else 'svd')
                if res < 1e-8: ok = True; break
            if ok: best = r; break
        out[q] = best
        print(('Bell' if bell else 'Full'), 'k', k, 'q', q, 'min rank found', best, flush=True)
    return out
if __name__ == '__main__':
    import sys
    k = int(sys.argv[1]); bell = sys.argv[2] == 'bell'; qs = [int(x) for x in sys.argv[3].split(',')]; starts = int(sys.argv[4])
    frontier(k, bell, 4**k, qs, starts)
