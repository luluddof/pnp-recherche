# Is SWAP (two k-qubit registers) Clifford-equivalent to M12 (CNOT layer)? search g in Sp(4k,2) with g F_SWAP g^-1 = F12
import numpy as np
from sympl_check import blocks, issymp, inv2, rand_symp
from sympl_check2 import nullspace2
for k in (1, 2, 3):
    F12, FM, Om = blocks(k); I = np.eye(k, dtype=np.int64); Z = np.zeros((k, k), dtype=np.int64)
    # SWAP of registers x <-> y: (ax, ay, bx, by) -> (ay, ax, by, bx)
    FS = np.block([[Z, I, Z, Z], [I, Z, Z, Z], [Z, Z, Z, I], [Z, Z, I, Z]])
    d = 4*k; rows = []
    for i in range(d):
        for j in range(d):
            r = np.zeros(d*d, dtype=np.int64)
            for l in range(d): r[i*d + l] ^= FS[l, j]
            for l in range(d): r[l*d + j] ^= F12[i, l]
            rows.append(r)
    basis = np.array(nullspace2(np.array(rows))); rng = np.random.default_rng(0); found = None
    for _ in range(200000):
        X = ((rng.integers(0, 2, len(basis)) @ basis) % 2).reshape(d, d)
        if issymp(X, Om): found = X; break
    ok = found is not None and np.array_equal((found @ FS @ inv2(found)) % 2, F12)
    print(f'k={k}: symplectic X with X F_SWAP X^-1 = F_12 found: {ok}; dim of linear solution space {len(basis)}')
