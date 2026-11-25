# rank of e_j(J_d) = sum_{|S|=j} J_d^{(x)S} (x) I^{(x)S^c} on (C^d)^{(x)m}  vs  sum_w min(dim_w, dim_{w-j})  (weight = sum of levels)
import numpy as np, itertools, math
from merged import rk, kron_list
def ej(d, m, j):
    J = np.zeros((d, d), dtype=np.int64)
    for a in range(d-1): J[a, a+1] = 1
    I = np.eye(d, dtype=np.int64); A = None
    for S in itertools.combinations(range(m), j):
        T = kron_list([J if s in S else I for s in range(m)]); A = T if A is None else A + T
    return A
def pred(d, m, j):
    dims = np.zeros(m*(d-1)+1, dtype=np.int64)
    for w in itertools.product(range(d), repeat=m): dims[sum(w)] += 1
    return sum(min(dims[w], dims[w-j]) for w in range(j, len(dims)))
for d, M in [(3, 6), (4, 5)]:
    for m in range(1, M+1):
        print(d, m, [(rk(ej(d, m, j))[0], pred(d, m, j)) for j in range(0, m+1)], flush=True)
