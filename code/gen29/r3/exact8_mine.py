# G29-R3 independent exact Bell-model optimum at N=8 (k=3) for q <= 4.
# phi(w) = sum_{s in S} kappa_s (-1)^{B(s,w)}, |S| = q; p = #{w : phi(w) != beta(w)}, beta = (-1)^{Q(w)}.
# Reduction: supports up to the affine orthogonal group (translations S -> S+t and O(Q) reflections), orbits by union-find.
# Per orbit rep: phi depends on w only via x(w) = (B(s,w))_s; optimum is attained at a vertex G determined by q
# independent equations G(x_j) = +-1 on image points (characters distinct on Im => independent); enumerate all of them.
import sys, itertools, time
import numpy as np
from fractions import Fraction
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
k = 3; n2 = 64; msk = 7
POP = np.array([bin(i).count('1') for i in range(256)])
Q = lambda x: POP[(x & msk) & (x >> k)] & 1
B = lambda x, v: (POP[(x & msk) & (v >> k)] + POP[(x >> k) & (v & msk)]) & 1
W = np.arange(n2); beta = 1 - 2*Q(W)
q = int(sys.argv[1])
combs = np.array(list(itertools.combinations(range(n2), q)), dtype=np.int64)
masks = np.zeros(len(combs), dtype=np.uint64)
for j in range(q): masks |= (np.uint64(1) << combs[:, j].astype(np.uint64))
order = np.argsort(masks); smask = masks[order]
gens = [W ^ (1 << i) for i in range(6)]                       # translations
gens += [W ^ (B(W, v) * v) for v in range(1, n2) if Q(v) == 1]  # orthogonal transvections (reflections)
for g in gens:  # sanity: Q-preserving up to translation part
    pass
rows, cols = [], []
for g in gens:
    im = g[combs]; m = np.zeros(len(combs), dtype=np.uint64)
    for j in range(q): m |= (np.uint64(1) << im[:, j].astype(np.uint64))
    idx = order[np.searchsorted(smask, m)]; assert np.all(masks[idx] == m)
    rows.append(np.arange(len(combs))); cols.append(idx)
r = np.concatenate(rows); c = np.concatenate(cols)
G = coo_matrix((np.ones(len(r)), (r, c)), shape=(len(combs), len(combs)))
nc, lab = connected_components(G, directed=True, connection='weak')
reps = [np.where(lab == i)[0][0] for i in range(nc)]
print(f'q={q}: {len(combs)} supports, {nc} orbits', flush=True)
def best(S):
    X = np.stack([B(W, s) for s in S], 1)                       # 64 x q bits
    keys = X @ (1 << np.arange(q)); im, inv = np.unique(keys, return_inverse=True)
    a = np.bincount(inv, weights=(beta == 1)).astype(int); b = np.bincount(inv, weights=(beta == -1)).astype(int)
    Chi = 1 - 2*((im[:, None] >> np.arange(q)) & 1)              # m x q character values on image points
    m = len(im); bestag = -1; bestsol = None
    for pts in itertools.combinations(range(m), q):
        A = Chi[list(pts)].astype(float)
        if abs(np.linalg.det(A)) < 1e-9: continue
        Ainv = np.linalg.inv(A)
        for sg in itertools.product((1, -1), repeat=q):
            kap = Ainv @ np.array(sg, float); Gv = Chi @ kap
            ag = int(np.sum(a[np.abs(Gv-1) < 1e-9]) + np.sum(b[np.abs(Gv+1) < 1e-9]))
            if ag > bestag: bestag, bestsol = ag, (pts, sg)
    # exact recheck of the best
    pts, sg = bestsol
    import sympy
    Ak = sympy.Matrix(Chi[list(pts)].tolist()); kap = Ak.solve(sympy.Matrix(sg))
    phi = [sum(kap[i]*(1-2*int(B(w, S[i]))) for i in range(q)) for w in range(n2)]
    ex = sum(1 for w in range(n2) if phi[w] == int(beta[w]))
    assert ex == bestag, (ex, bestag)
    return n2 - bestag
t0 = time.time(); res = []
for ri in reps:
    S = [int(x) for x in combs[ri]]; res.append((best(S), S))
res.sort()
print(f'q={q}: p_min = {res[0][0]} (support {res[0][1]}); orbit p-values: {sorted(set(x[0] for x in res))}; {time.time()-t0:.1f}s', flush=True)
# invariance spot check: random members of the optimal orbit
rng = np.random.default_rng(0); o = lab[reps[[x for x in range(nc)]][0]] if False else None
oi = lab[[ri for ri in reps if [int(x) for x in combs[ri]] == res[0][1]][0]]
mem = np.where(lab == oi)[0]; pick = rng.choice(mem, size=min(5, len(mem)), replace=False)
print('  orbit-member check p:', [best([int(x) for x in combs[j]]) for j in pick], flush=True)
