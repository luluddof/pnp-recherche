# numeric check: Pauli-sparse approximation of SWAP_N (N=4, 2+2 qubits) gives the same frontier as for CZ^{(x)2}
import numpy as np
from core import *
from sl import sl
k = 2; N = 4; D = 16
SW = np.zeros((D, D))
for x in range(N):
    for y in range(N): SW[y + N*x, x + N*y] = 1     # index = low k bits register 1, high bits register 2
rng = np.random.default_rng(4); mask = np.ones((D, D))
for q in (1, 2, 3, 4, 5, 6):
    for r in range(int(np.ceil(D/(4*q))), D+1):
        ok = any(sl(SW, q, r, mask, rng, 250, init='rand' if s else 'svd')[0] < 1e-8 for s in range(8))
        if ok: print('SWAP N=4 Pauli-sparse: q', q, 'min rank found', r, flush=True); break
