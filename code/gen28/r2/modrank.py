import numpy as np
def rank_mod(M, p):
    A = np.array(M, dtype=np.int64) % p
    n, m = A.shape; r = 0
    for c in range(m):
        if r == n: break
        piv = np.nonzero(A[r:, c])[0]
        if len(piv) == 0: continue
        i = r + piv[0]
        if i != r: A[[r, i]] = A[[i, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        rows = np.nonzero(A[:, c])[0]; rows = rows[rows != r]
        if len(rows):
            f = A[rows, c].copy()
            A[rows] = (A[rows] - (f[:, None] * A[r][None, :]) % p) % p
        r += 1
    return r
