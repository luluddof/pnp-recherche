# core tools: real Paulis X^a Z^b on n qubits (bitmask a,b), Pauli transform, target D_beta = CZ^{(x)k}
import numpy as np
def popc(x): return bin(x).count('1')
def fwht(a, axis=-1):
    a = np.array(a, dtype=float, copy=True); a = np.moveaxis(a, axis, -1); n = a.shape[-1]; h = 1
    while h < n:
        a = a.reshape(a.shape[:-1] + (n//(2*h), 2, h))
        x = a[..., 0, :].copy(); y = a[..., 1, :].copy()
        a[..., 0, :] = x + y; a[..., 1, :] = x - y
        a = a.reshape(a.shape[:-3] + (n,)); h *= 2
    return np.moveaxis(a, -1, axis)
def beta(k):
    # z = (x bits 0..k-1, y bits k..2k-1); beta(z) = (-1)^{x.y}  -> D_beta = CZ^{(x)k} on pairs (i, k+i)
    D = 4**k; z = np.arange(D); x = z & ((1<<k)-1); y = z >> k
    return np.array([(-1)**popc(int(a & b)) for a, b in zip(x, y)], dtype=float)
def pauli(a, b, n):
    D = 1<<n; z = np.arange(D); M = np.zeros((D, D))
    sgn = np.array([(-1)**popc(int(b & zz)) for zz in z])
    M[z ^ a, z] = sgn            # X^a Z^b |z> = (-1)^{b.z} |z+a>
    return M
def ptrans(M):
    # coefficients c[a,b] with M = sum c[a,b] X^a Z^b ; d_a(z) = M[z^a, z] = sum_b c[a,b] (-1)^{b.z}
    D = M.shape[0]; z = np.arange(D)
    d = np.array([M[z ^ a, z] for a in range(D)])      # D x D
    return fwht(d, axis=1) / D
def pinv_trans(c):
    D = c.shape[0]; z = np.arange(D); d = fwht(c, axis=1)  # d_a(z)
    M = np.zeros((D, D))
    for a in range(D): M[z ^ a, z] = d[a]
    return M
def nrank(M, tol=1e-8):
    s = np.linalg.svd(M, compute_uv=False); return int((s > tol * max(1, s[0])).sum()), s
