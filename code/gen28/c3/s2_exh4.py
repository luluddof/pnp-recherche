# exhaustive k<=4: max mu_{1/4} agreement of degree<=d multilinear g with parity
import numpy as np, itertools
from math import comb
for k in [3, 4]:
    pats = np.array(list(itertools.product([0,1], repeat=k)))
    mu = np.prod(np.where(pats == 1, 0.25, 0.75), axis=1); par = (-1.0)**pats.sum(1)
    P = len(pats)
    for d in range(k+1):
        mons = [U for r in range(d+1) for U in itertools.combinations(range(k), r)]
        A = np.array([[np.prod(T[list(U)]) if U else 1.0 for U in mons] for T in pats])
        best = 0; bestZ = None
        for mask in range(1, 2**P):
            Z = [i for i in range(P) if mask >> i & 1]
            m = mu[Z].sum()
            if m <= best + 1e-12: continue
            AZ = A[Z]; r1 = np.linalg.matrix_rank(AZ); r2 = np.linalg.matrix_rank(np.hstack([AZ, par[Z, None]]))
            if r1 == r2: best = m; bestZ = Z
        symw = None
        print(f"k={k} d={d} exhaustive max agreement = {best:.4f}  Z={[''.join(map(str,pats[i])) for i in bestZ]}", flush=True)
