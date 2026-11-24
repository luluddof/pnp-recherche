# Highest-weight-vector interpolation.  For a weight lam = (l1,l2,l3,l4) (partitions, <= 3 rows) the h.w. vectors of weight lam in
# S^d(W*) are spanned by "tableau polynomials" f_T (column antisymmetrisations, one filling per leg; leg 1 canonical).
# We pick a basis (rank on generic points must reach g(lam) = generalised Kronecker coefficient), evaluate on random points of
# J(p,q) and of Z3 (mod PR).  eqs = g - rank_J ; sep = rank_{J+Z3} - rank_J  (> 0 <=> some equation of J of this weight is nonzero
# on Z3, i.e. separates the swap).
import numpy as np, itertools, sys, time
from common import *
from kron import g as kron_g
def cells(lam):
    return [(r, j) for r, L in enumerate(lam) for j in range(L)]
def leg_assignments(lam, filling):
    """filling: list cell index for copy c -> returns array (n_assign, d) of indices and signs"""
    cl = cells(lam); d = len(cl)
    ncol = lam[0]; cols = [[k for k, (r, j) in enumerate(cl) if j == jj] for jj in range(ncol)]
    cell_of_copy = filling   # copy c sits in cell filling[c]
    copy_of_cell = {cell_of_copy[c]: c for c in range(d)}
    per_col = []
    for col in cols:
        h = len(col); opts = []
        for perm in itertools.permutations(range(h)):
            sgn = 1
            for i in range(h):
                for j in range(i+1, h):
                    if perm[i] > perm[j]: sgn = -sgn
            opts.append((perm, sgn))
        per_col.append((col, opts))
    out_idx = []; out_sg = []
    for combo in itertools.product(*[o for _, o in per_col]):
        idx = [0]*d; s = 1
        for (col, _), (perm, sg) in zip(per_col, combo):
            for k, cellk in enumerate(col):
                r = cl[cellk][0]; idx[copy_of_cell[cellk]] = perm[r]
            s *= sg
        out_idx.append(idx); out_sg.append(s)
    return np.array(out_idx, dtype=np.int64), np.array(out_sg, dtype=np.int64)
def terms(lams, fills, N=3):
    legs = [leg_assignments(l, f) for l, f in zip(lams, fills)]
    var = np.zeros((1, len(fills[0])), dtype=np.int64); sg = np.ones(1, dtype=np.int64)
    for (idx, s), w in zip(legs, [27, 9, 3, 1]):
        var = (var[:, None, :] + w*idx[None, :, :]).reshape(-1, idx.shape[1]); sg = (sg[:, None]*s[None, :]).ravel()
    return var, sg
def evalT(var, sg, P, chunk=200000):
    """P: (npts, 81) mod PR -> values (npts,)"""
    tot = np.zeros(P.shape[0], dtype=np.int64)
    for s in range(0, len(sg), chunk):
        v = var[s:s+chunk]; prod = np.ones((P.shape[0], len(v)), dtype=np.int64)
        for c in range(v.shape[1]): prod = (prod * P[:, v[:, c]]) % PR
        tot = (tot + (prod * (sg[s:s+chunk] % PR)[None, :]) % PR @ np.ones(len(v), dtype=np.int64)) % PR
    return tot
def analyse(lams, p=2, q=2, seed=0, budget=3, extra=25, verbose=True):
    lams = [tuple(x for x in l if x > 0) for l in lams]
    d = sum(lams[0]); gl = kron_g(*lams); rng = np.random.default_rng(seed)
    if gl == 0: return dict(g=0)
    ntest = gl + extra
    Pgen = np.array([rnd(81, rng) for _ in range(ntest)]); PJ = np.array([ptJ(3,p,q,rng).ravel() for _ in range(ntest)])
    PZ = np.array([ptZ3(3,rng).ravel() for _ in range(20)])
    can = list(range(d)); basis = []; Mg = np.zeros((0, ntest), dtype=np.int64); r = 0; tries = 0; t0 = time.time()
    while r < gl and tries < budget*gl + 20:
        tries += 1
        fills = [can] + [list(rng.permutation(d)) for _ in range(3)]
        var, sg = terms(lams, fills)
        vg = evalT(var, sg, Pgen)
        M2 = np.vstack([Mg, vg[None, :]])
        r2 = rank(M2)
        if r2 > r: Mg = M2; r = r2; basis.append(fills)
    if r < gl: return dict(g=gl, found=r, error='basis incomplete')
    EJ = []; EZ = []
    for fills in basis:
        var, sg = terms(lams, fills); EJ.append(evalT(var, sg, PJ)); EZ.append(evalT(var, sg, PZ))
    EJ = np.array(EJ).T; EZ = np.array(EZ).T      # rows: points, cols: basis polys
    rJ = rank(EJ); rJZ = rank(np.vstack([EJ, EZ]))
    return dict(g=gl, eqs=gl - rJ, sep=rJZ - rJ, nterms=len(sg), basis=basis, time=time.time()-t0, EJ=EJ, EZ=EZ)
if __name__ == "__main__":
    import json
    spec = json.loads(sys.argv[1])   # list of [mu, nu] pairs; weight (mu,nu,nu,mu)
    p, q = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (2, 2)
    for mu, nu in spec:
        res = analyse([mu, nu, nu, mu], p, q)
        print(tuple(mu), tuple(nu), {k: v for k, v in res.items() if k not in ('basis','EJ','EZ')}, 'WITNESS' if res.get('sep', 0) > 0 else '', flush=True)
        if res.get('sep', 0) > 0:
            np.save(f"witness_{''.join(map(str,mu))}_{''.join(map(str,nu))}.npy", np.array(res['basis'], dtype=object), allow_pickle=True)
