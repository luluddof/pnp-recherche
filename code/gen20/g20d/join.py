# Terracini: dim of the join J(p,q) = closure{A + B : rank_1 A <= p, rank_2 B <= q} in End(K^N (x) K^N), exact mod a prime.
import numpy as np, sys, random
from flint import nmod_mat
PR = 1000003; random.seed(1)
def rank(rows, ncol):
    return nmod_mat(len(rows), ncol, [v % PR for r in rows for v in r], PR).rank() if rows else 0
def tangent_rows(N, p, s):
    D = N*N
    U = [[random.randrange(PR) for _ in range(D)] for _ in range(p)]
    W = [[random.randrange(PR) for _ in range(D)] for _ in range(p)]
    rows = []
    for i in range(p):
        for w in range(D):
            M = np.zeros((D,D), dtype=object); M[:, w] = U[i]; rows.append(M)      # U_i e_w^T
            M = np.zeros((D,D), dtype=object); M[w, :] = W[i]; rows.append(M)      # e_w W_i^T
    out = []
    for M in rows:
        T = M.reshape(N,N,N,N)            # indices of the flattening s
        if s == 2: T = T.transpose(0,2,1,3)  # (x,u,y,v) -> (x,y,u,v)
        out.append([int(v) for v in T.flatten()])
    return out
for N in [2,3,4]:
    D = N*N
    res = {}
    for p in range(0, D+1):
        Rp = tangent_rows(N, p, 1)
        for q in range(0, D+1):
            if p + q == 0: continue
            if (p, q) in res: continue
            dim = rank(Rp + tangent_rows(N, q, 2), D*D)
            exp = min(D*D, (2*p*D - p*p) + (2*q*D - q*q))
            res[(p,q)] = (dim, exp)
            if dim == D*D: break
    fills = sorted([(p,q) for (p,q),(d,e) in res.items() if d == D*D])
    front = [(p, min(q for (pp,q) in fills if pp == p)) for p in sorted(set(p for p,_ in fills))]
    defect = [(p,q,d,e) for (p,q),(d,e) in res.items() if d != e]
    print('N',N,'generic frontier (least q filling, per p):', front, flush=True)
    print('   defective (p,q,dim,expected):', defect[:20], flush=True)
