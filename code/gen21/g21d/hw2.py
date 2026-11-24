# Faster evaluation of tableau polynomials: sum over legs 1-3 column-permutations explicitly, leg 4 via determinants.
import numpy as np, sys, time, json
from common import *
from kron import g as kron_g
from hw import cells, leg_assignments
def leg4_cols(lam, filling):
    cl = cells(lam); copy_of_cell = {filling[c]: c for c in range(len(filling))}
    cols = []
    for jj in range(lam[0]):
        cc = sorted([k for k, (r, j) in enumerate(cl) if j == jj], key=lambda k: cl[k][0])
        cols.append([copy_of_cell[k] for k in cc])
    return cols
def det_mod(Ms, h):   # Ms: (..., h, h) int64 mod PR
    if h == 1: return Ms[..., 0, 0] % PR
    if h == 2: return (Ms[...,0,0]*Ms[...,1,1] - Ms[...,0,1]*Ms[...,1,0]) % PR
    a = Ms
    t1 = a[...,0,0] * ((a[...,1,1]*a[...,2,2] - a[...,1,2]*a[...,2,1]) % PR) % PR
    t2 = a[...,0,1] * ((a[...,1,0]*a[...,2,2] - a[...,1,2]*a[...,2,0]) % PR) % PR
    t3 = a[...,0,2] * ((a[...,1,0]*a[...,2,1] - a[...,1,1]*a[...,2,0]) % PR) % PR
    return (t1 - t2 + t3) % PR
def prep(lams, fills):
    legs = [leg_assignments(l, f) for l, f in zip(lams[:3], fills[:3])]
    d = len(fills[0]); var = np.zeros((1, d), dtype=np.int64); sg = np.ones(1, dtype=np.int64)
    for (idx, s), w in zip(legs, [9, 3, 1]):
        var = (var[:, None, :] + w*idx[None, :, :]).reshape(-1, d); sg = (sg[:, None]*s[None, :]).ravel()
    return var, sg % PR, leg4_cols(lams[3], fills[3])
def evalT2(pre, P, chunk=20000):
    var, sg, cols = pre; npts = P.shape[0]; P3 = P.reshape(npts, 27, 3); tot = np.zeros(npts, dtype=np.int64)
    for s in range(0, len(sg), chunk):
        v = var[s:s+chunk]; acc = np.ones((npts, len(v)), dtype=np.int64)
        for col in cols:
            h = len(col)
            M = np.stack([P3[:, v[:, c], :h] for c in col], axis=2)   # (npts, nv, h rows, h cols)
            acc = (acc * det_mod(M, h)) % PR
        tot = (tot + (acc * sg[s:s+chunk][None, :] % PR).sum(axis=1)) % PR
    return tot
def analyse(lams, p=2, q=2, seed=0, budget=3, extra=25):
    lams = [tuple(x for x in l if x > 0) for l in lams]
    d = sum(lams[0]); gl = kron_g(*lams); rng = np.random.default_rng(seed)
    if gl == 0: return dict(g=0)
    ntest = gl + extra
    Pgen = np.array([rnd(81, rng) for _ in range(ntest)]); PJ = np.array([ptJ(3,p,q,rng).ravel() for _ in range(ntest)])
    PZ = np.array([ptZ3(3,rng).ravel() for _ in range(20)])
    can = list(range(d)); basis = []; vals = []; r = 0; tries = 0; t0 = time.time()
    while r < gl and tries < budget*gl + 20:
        tries += 1
        fills = [can] + [list(rng.permutation(d)) for _ in range(3)]
        pre = prep(lams, fills); vg = evalT2(pre, Pgen)
        r2 = rank(np.array(vals + [vg]))
        if r2 > r: vals.append(vg); r = r2; basis.append(fills)
    if r < gl: return dict(g=gl, found=r, error='basis incomplete')
    EJ = []; EZ = []
    for fills in basis:
        pre = prep(lams, fills); EJ.append(evalT2(pre, PJ)); EZ.append(evalT2(pre, PZ))
    EJ = np.array(EJ).T; EZ = np.array(EZ).T
    rJ = rank(EJ); rJZ = rank(np.vstack([EJ, EZ]))
    return dict(g=gl, eqs=gl - rJ, sep=rJZ - rJ, nterms=len(pre[1]), tries=tries, time=round(time.time()-t0,1), basis=basis, EJ=EJ, EZ=EZ)
if __name__ == "__main__":
    spec = json.loads(sys.argv[1]); p, q = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (2, 2)
    for lam in spec:
        if len(lam) == 2: lam = [lam[0], lam[1], lam[1], lam[0]]
        res = analyse(lam, p, q)
        print([tuple(l) for l in lam], {k: v for k, v in res.items() if k not in ('basis','EJ','EZ')}, 'WITNESS' if res.get('sep', 0) > 0 else '', flush=True)
        if res.get('eqs', 0) > 0:
            np.save('eqs_' + '_'.join(''.join(map(str,l)) for l in lam) + '.npy', np.array(res['basis'], dtype=object), allow_pickle=True)
