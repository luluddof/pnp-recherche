"""Scaling of poly-time-findable short even covers (H1 cycle search) at m = c*n^1.4 vs c*n^1.5,
and the needed LP value B12/2 (computed with sparse eigsh). Upper bound on any cover-LP: (#clauses in covers)/minsize."""
import sys, random, math, numpy as np, scipy.sparse as sp, time
from scipy.sparse.linalg import eigsh
from fko import gen, cycle_covers, check_even
def B12(n, V, S):
    c = np.zeros(n); r, co, va = [], [], []
    for v, s in zip(V, S):
        for i in range(3): c[v[i]] += s[i]
        for i, j in ((0,1),(0,2),(1,2)):
            r += [v[i], v[j]]; co += [v[j], v[i]]; va += [-s[i]*s[j]/2]*2
    A = sp.csr_matrix((va, (r, co)), shape=(n, n))
    return np.abs(c).sum() + n*max(0, eigsh(A, k=1, which='LA', return_eigenvectors=False)[0])
expo = float(sys.argv[1]); c = float(sys.argv[2])
for n in [int(x) for x in sys.argv[3].split(',')]:
    m = int(c*n**expo); rng = random.Random(n); V, S = gen(n, m, rng)
    t = time.time(); covs, core = cycle_covers(n, V, 4, rng); tc = time.time()-t
    sgn = np.prod(S, axis=1); bad = [K for K in covs if np.prod([sgn[k] for k in K]) == -1]
    used = set().union(*bad) if bad else set()
    sz = sorted(len(K) for K in bad)
    need = B12(n, V, S)/2
    print(f"n={n} m={m} m/n^1.5={m/n**1.5:.2f} core={core} nbad={len(bad)} minsz={sz[0] if sz else None} medsz={sz[len(sz)//2] if sz else None} clauses_in_badcovers={len(used)} need={need:.0f} trivialLPmax~{len(used)/max(1,(sz[0] if sz else 1)):.0f} t={tc:.1f}s", flush=True)
