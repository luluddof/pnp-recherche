# Cayley-Dickson sign tables: e_i e_j = s(i,j) e_{i^j}, dim N = 2^n.  (a,b)(c,d) = (ac - conj(d) b, d a + b conj(c)).
from functools import lru_cache
import numpy as np

@lru_cache(maxsize=None)
def mul(i, j, n):
    if n == 0: return 1, 0
    h = 1 << (n - 1)
    if i < h and j < h:
        s, k = mul(i, j, n - 1); return s, k
    if i < h and j >= h:          # (e_i,0)(0,e_j') = (0, e_j' e_i)
        s, k = mul(j - h, i, n - 1); return s, k + h
    if i >= h and j < h:          # (0,e_i')(e_j,0) = (0, e_i' conj(e_j))
        s, k = mul(i - h, j, n - 1); c = 1 if j == 0 else -1; return s * c, k + h
    # (0,e_i')(0,e_j') = (-conj(e_j') e_i', 0)
    s, k = mul(j - h, i - h, n - 1); c = 1 if j == h else -1; return -s * c, k

def table(n):
    N = 1 << n
    S = np.zeros((N, N), dtype=np.int64); K = np.zeros((N, N), dtype=np.int64)
    for i in range(N):
        for j in range(N):
            S[i, j], K[i, j] = mul(i, j, n)
    return S, K

if __name__ == '__main__':
    for n in range(5):
        S, K = table(n); N = 1 << n
        assert all((K[i, j] == (i ^ j)) for i in range(N) for j in range(N))
        # norm multiplicativity check on random integer vectors
        rng = np.random.default_rng(1); bad = 0
        for t in range(20):
            a = rng.integers(-3, 4, N); b = rng.integers(-3, 4, N); c = np.zeros(N, dtype=np.int64)
            for i in range(N):
                for j in range(N): c[K[i, j]] += S[i, j] * a[i] * b[j]
            bad += (c @ c != (a @ a) * (b @ b))
        print(n, N, 'composition failures', bad, '/20')
