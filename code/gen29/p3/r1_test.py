# R_V(1) = n - max eigenvalue multiplicity of P^dagger V, over all 2k-qubit Paulis P; V = M12 (CNOT layer) and M.
import numpy as np, sys, itertools
from core import Mperm, pauli_M
k = int(sys.argv[1]); n = 4**k; N = 2**k; msk = N-1
M = Mperm(k).astype(float)
M12 = np.zeros((n, n))
for i in range(n):
    x, y = i & msk, i >> k; M12[(x ^ y) | (y << k), i] = 1
def maxmult(A):
    ev = np.linalg.eigvals(A); best = 0
    for e in ev:
        best = max(best, int(np.sum(np.abs(ev - e) < 1e-6)))
    return best
for name, V in (('M12', M12), ('M', M)):
    res = {}
    for a1, a2, c, d in itertools.product(range(N), repeat=4):
        P = pauli_M(k, a1, a2, c, d).astype(float)
        inL12 = (a2 == 0 and d == 0)
        m = maxmult(P.T @ V); key = 'L12' if inL12 else 'other'
        res[key] = max(res.get(key, 0), m)
    print(f'k={k} {name}: R(1) over L12 Paulis = {n - res["L12"]}, over other Paulis = {n - res["other"]}')
