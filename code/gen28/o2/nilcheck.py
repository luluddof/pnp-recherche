# check the inclusion-matrix rank formula for Jordan 2-blocks: M = X (J_2^{+p} + 0) X^-1, class i = #C factors (I-M)
import numpy as np, itertools, math
from merged import rk, kron_list
def cls(M, k, i):
    n = M.shape[0]; I = np.eye(n, dtype=np.int64); A = None
    for Q in itertools.combinations(range(k), i):
        T = kron_list([(I - M) if s in Q else M for s in range(k)]); A = T if A is None else A + T
    return A
def formula(n2, p, k, j):
    t = 0
    for m in range(j, k+1):
        inc = sum(min(math.comb(m, w), math.comb(m, w-j)) if w-j >= 0 else 0 for w in range(m+1))
        t += math.comb(k, m)*p**m*(n2-2*p)**(k-m)*inc
    return t
rng = np.random.default_rng(0)
for n2, p in [(4, 1), (4, 2), (5, 2)]:
    J = np.zeros((n2, n2), dtype=np.int64)
    for b in range(p): J[2*b, 2*b+1] = 1
    X = np.eye(n2, dtype=np.int64) + np.triu(rng.integers(-2, 3, (n2, n2)), 1); Xi = np.rint(np.linalg.inv(X)).astype(np.int64)
    M = X @ J @ Xi
    for k in (2, 3):
        print(n2, p, k, [(rk(cls(M, k, i))[0], formula(n2, p, k, k-i)) for i in range(k+1)])
