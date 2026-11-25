# Exhaustive Bell model at N=4: is there phi in span{chi_s : s in S}, |S|=q, agreeing with beta on a set Z of size z?
import numpy as np, itertools, sys
pts = np.array(list(itertools.product([0, 1], repeat=4)))  # v=(c1,c2,d1,d2)
beta = (-1.0)**((pts[:, 0]*pts[:, 2] + pts[:, 1]*pts[:, 3]) % 2)
# characters chi_(a,b)(c,d) = (-1)^{a.d+b.c}
lab = pts
CH = (-1.0)**((lab[:, :2] @ pts[:, 2:].T + lab[:, 2:] @ pts[:, :2].T) % 2)   # [label, point]
def feasible(q, z):
    Zs = np.array(list(itertools.combinations(range(16), z)))
    for S in itertools.combinations(range(16), q):
        A = CH[list(S)].T            # 16 x q
        AZ = A[Zs]                   # nZ x z x q
        bZ = beta[Zs][:, :, None]
        r1 = np.linalg.matrix_rank(AZ); r2 = np.linalg.matrix_rank(np.concatenate([AZ, bZ], 2))
        hit = np.flatnonzero(r1 == r2)
        if hit.size: return S, Zs[hit[0]]
    return None
for q, p in [(1, 6), (3, 6), (4, 4), (5, 4), (6, 1), (15, 1)]:
    up = feasible(q, 16 - p)          # achievable p
    lo = feasible(q, 16 - p + 1)      # p-1 achievable?
    print(f"q={q}: p={p} achievable: {up is not None}; p={p-1} achievable: {lo is not None}", flush=True)
