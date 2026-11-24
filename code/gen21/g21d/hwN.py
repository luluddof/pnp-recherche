# hw2.py for general N (N = 4 tests).  Same method: tableau h.w. vectors, legs 1-3 explicit, leg 4 by determinants (h <= 4).
import numpy as np, sys, time, json
from common import *
from kron import g as kron_g
from hw import cells, leg_assignments
from hw2 import leg4_cols
def det_mod(M, h):
    if h == 1: return M[..., 0, 0] % PR
    tot = 0
    for j in range(h):   # Laplace along first row
        minor = np.delete(np.delete(M, 0, axis=-2), j, axis=-1)
        tot = (tot + (-1)**j * (M[..., 0, j] * det_mod(minor, h-1) % PR)) % PR
    return tot % PR
def prep(lams, fills, N):
    legs = [leg_assignments(l, f) for l, f in zip(lams[:3], fills[:3])]
    d = len(fills[0]); var = np.zeros((1, d), dtype=np.int64); sg = np.ones(1, dtype=np.int64)
    for (idx, s), w in zip(legs, [N*N, N, 1]):
        var = (var[:, None, :] + w*idx[None, :, :]).reshape(-1, d); sg = (sg[:, None]*s[None, :]).ravel()
    return var, sg % PR, leg4_cols(lams[3], fills[3])
def evalT(pre, P, N, chunk=5000):
    var, sg, cols = pre; npts = P.shape[0]; P3 = P.reshape(npts, N**3, N); tot = np.zeros(npts, dtype=np.int64)
    for s in range(0, len(sg), chunk):
        v = var[s:s+chunk]; acc = np.ones((npts, len(v)), dtype=np.int64)
        for col in cols:
            h = len(col); M = np.stack([P3[:, v[:, c], :h] for c in col], axis=2)
            acc = (acc * det_mod(M, h)) % PR
        tot = (tot + (acc * sg[s:s+chunk][None, :] % PR).sum(axis=1)) % PR
    return tot
def analyse(lams, N, p, q, seed=0, budget=3, extra=25):
    lams = [tuple(x for x in l if x > 0) for l in lams]
    d = sum(lams[0]); gl = kron_g(*lams); rng = np.random.default_rng(seed)
    if gl == 0: return dict(g=0)
    ntest = gl + extra
    Pgen = np.array([rnd(N**4, rng) for _ in range(ntest)]); PJ = np.array([ptJ(N,p,q,rng).ravel() for _ in range(ntest)])
    PZ = np.array([ptZ3(N,rng).ravel() for _ in range(20)])
    can = list(range(d)); basis = []; vals = []; r = 0; tries = 0; t0 = time.time()
    while r < gl and tries < budget*gl + 20:
        tries += 1
        fills = [can] + [list(rng.permutation(d)) for _ in range(3)]
        vg = evalT(prep(lams, fills, N), Pgen, N)
        r2 = rank(np.array(vals + [vg]))
        if r2 > r: vals.append(vg); r = r2; basis.append(fills)
    if r < gl: return dict(g=gl, found=r, error='basis incomplete')
    EJ = np.array([evalT(prep(lams, f, N), PJ, N) for f in basis]).T; EZ = np.array([evalT(prep(lams, f, N), PZ, N) for f in basis]).T
    rJ = rank(EJ); rJZ = rank(np.vstack([EJ, EZ]))
    return dict(g=gl, eqs=gl - rJ, sep=rJZ - rJ, tries=tries, time=round(time.time()-t0,1))
if __name__ == "__main__":
    N, p, q = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    for lam in json.loads(sys.argv[1]):
        if len(lam) == 2: lam = [lam[0], lam[1], lam[1], lam[0]]
        res = analyse(lam, N, p, q)
        print(N, (p,q), [tuple(l) for l in lam], res, 'WITNESS' if res.get('sep', 0) > 0 else '', flush=True)
