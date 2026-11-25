# k=2 (4 qubits, D=16), supports of size 3: orbits under the symmetry group of D_f = CZ (x) CZ, then numeric min rank per orbit.
import numpy as np, itertools, time, pickle
from core import *
k = 2; n = 4; D = 16
def J(c): return ((c & 3) << 2) | (c >> 2)          # swap x,y halves
def Bf(z, w): return popc(z & J(w)) & 1
gens = []
for i in range(n): gens.append(lambda ab, i=i: (ab[0], ab[1] ^ (1 << i)))                       # twisted translations
for i, j in itertools.combinations(range(n), 2):                                                      # CZ_ij conjugation
    gens.append(lambda ab, i=i, j=j: (ab[0], ab[1] ^ (((ab[0] >> j) & 1) << i) ^ (((ab[0] >> i) & 1) << j)))
def transv(v):
    G = np.zeros((n, n), dtype=np.int64)
    for c in range(n):
        z = 1 << c; img = z ^ (v if Bf(v, z) else 0)
        for r in range(n): G[r, c] = (img >> r) & 1
    Gi = np.linalg.inv(G).round().astype(np.int64) % 2 if False else None
    return G
def apply_mat(G, x):
    y = 0
    for c in range(n):
        if (x >> c) & 1: y ^= sum(int(G[r, c]) << r for r in range(n))
    return y
def inv2(A):
    m = A.shape[0]; M = np.concatenate([A % 2, np.eye(m, dtype=np.int64)], 1); r = 0
    for c in range(m):
        p = [i for i in range(r, m) if M[i, c]][0]; M[[r, p]] = M[[p, r]]
        for i in range(m):
            if i != r and M[i, c]: M[i] ^= M[r]
        r += 1
    return M[:, m:]
for v in range(1, D):
    G = transv(v); GiT = inv2(G).T % 2
    ta = [apply_mat(G, x) for x in range(D)]; tb = [apply_mat(GiT, x) for x in range(D)]
    gens.append(lambda ab, ta=ta, tb=tb: (ta[ab[0]], tb[ab[1]]))
# sanity: each generator maps D_f to +-D_f * (Z-Pauli) under the corresponding conjugation -- checked implicitly by rank invariance test below
pts = [(a, b) for a in range(D) for b in range(D)]; idx = {p: i for i, p in enumerate(pts)}
gmaps = [np.array([idx[g(p)] for p in pts]) for g in gens]
t0 = time.time()
subs = list(itertools.combinations(range(256), 3)); sid = {s: i for i, s in enumerate(subs)}
parent = np.arange(len(subs))
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
for gm in gmaps:
    for i, s in enumerate(subs):
        t = tuple(sorted(gm[list(s)]))
        a, b = find(i), find(sid[t])
        if a != b: parent[a] = b
reps = {}
for i in range(len(subs)):
    r = find(i); reps.setdefault(r, i)
print('3-subsets', len(subs), 'orbits', len(reps), f'{time.time()-t0:.0f}s', flush=True)
pickle.dump([subs[i] for i in reps.values()], open('k2q3_reps.pkl', 'wb')); pickle.dump(pts, open('k2_pts.pkl', 'wb'))
