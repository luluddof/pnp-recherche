import numpy as np
def rank_mod(A, p):
    A = np.array(A, dtype=np.int64) % p
    A = A.copy(); n, m = A.shape; r = 0
    for c in range(m):
        if r == n: break
        piv = np.nonzero(A[r:, c])[0]
        if len(piv) == 0: continue
        i = r + piv[0]
        if i != r: A[[r, i]] = A[[i, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        col = A[:, c].copy(); col[r] = 0
        nz = np.nonzero(col)[0]
        if len(nz):
            A[nz] = (A[nz] - np.outer(col[nz], A[r])) % p
        r += 1
    return r
def solve_mod(M, b, p):
    """return one solution x of M x = b mod p, or None"""
    M = np.array(M, dtype=np.int64) % p; b = np.array(b, dtype=np.int64) % p
    n, m = M.shape
    A = np.concatenate([M, b[:, None]], axis=1)
    r = 0; pivs = []
    for c in range(m):
        if r == n: break
        piv = np.nonzero(A[r:, c])[0]
        if len(piv) == 0: continue
        i = r + piv[0]
        if i != r: A[[r, i]] = A[[i, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        col = A[:, c].copy(); col[r] = 0
        nz = np.nonzero(col)[0]
        if len(nz): A[nz] = (A[nz] - np.outer(col[nz], A[r])) % p
        pivs.append(c); r += 1
    if np.any(A[r:, m] % p): return None
    x = np.zeros(m, dtype=np.int64)
    for i, c in enumerate(pivs): x[c] = A[i, m]
    return x
