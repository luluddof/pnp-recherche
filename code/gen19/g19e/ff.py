# exact linear algebra mod p (batched), written for G19e
import numpy as np
def brank(A, p):
    """A: (B,n,m) int array; returns ranks mod p, batched Gaussian elimination."""
    A = (A % p).astype(np.int64).copy()
    B, n, m = A.shape
    rank = np.zeros(B, dtype=np.int64)
    idx = np.arange(B)
    inv = np.array([0] + [pow(a, p - 2, p) for a in range(1, p)], dtype=np.int64)
    for c in range(m):
        # find pivot row >= rank with nonzero in column c
        col = A[:, :, c].copy()
        mask = np.arange(n)[None, :] >= rank[:, None]
        nz = (col != 0) & mask
        has = nz.any(axis=1)
        piv = np.argmax(nz, axis=1)
        if not has.any():
            continue
        h = idx[has]
        r = np.minimum(rank[h], n - 1)
        # swap rows r and piv
        tmp = A[h, r, :].copy()
        A[h, r, :] = A[h, piv[h], :]
        A[h, piv[h], :] = tmp
        pr = A[h, r, :]
        pr = (pr * inv[pr[:, c]][:, None]) % p
        A[h, r, :] = pr
        f = A[h, :, c].copy()
        f[np.arange(len(h)), r] = 0
        A[h] = (A[h] - f[:, :, None] * pr[:, None, :]) % p
        rank[h] += 1
    return rank
def rank1(A, p):
    return int(brank(np.asarray(A)[None], p)[0])
def nullspace(A, p):
    """basis (rows) of the right kernel of A mod p"""
    A = (np.asarray(A) % p).astype(np.int64).copy()
    n, m = A.shape
    piv = []
    r = 0
    for c in range(m):
        nzr = [i for i in range(r, n) if A[i, c] % p]
        if not nzr: continue
        i = nzr[0]
        A[[r, i]] = A[[i, r]]
        A[r] = (A[r] * pow(int(A[r, c]), p - 2, p)) % p
        for j in range(n):
            if j != r and A[j, c]:
                A[j] = (A[j] - A[j, c] * A[r]) % p
        piv.append(c); r += 1
        if r == n: break
    free = [c for c in range(m) if c not in piv]
    basis = []
    for f in free:
        v = np.zeros(m, dtype=np.int64); v[f] = 1
        for i, c in enumerate(piv):
            v[c] = (-A[i, f]) % p
        basis.append(v)
    return np.array(basis, dtype=np.int64).reshape(len(basis), m)
