# own check of P3's group fact: F_12 = (g F_M g^-1)(h F_M h^-1) in Sp(4k,2), k = 2 (search), k = 4 (direct sum of k=2)
import numpy as np, sys, time
def blocks(k):
    I = np.eye(k, dtype=np.int64); Z = np.zeros((k, k), dtype=np.int64)
    # coordinates (ax, ay, bx, by), column vectors; F acts on columns
    F12 = np.block([[I, I, Z, Z], [Z, I, Z, Z], [Z, Z, I, Z], [Z, Z, I, I]]) % 2
    FM = np.block([[Z, I, Z, Z], [I, I, Z, Z], [Z, Z, I, I], [Z, Z, I, Z]]) % 2
    n = 2*k; Om = np.block([[np.zeros((n, n), dtype=np.int64), np.eye(n, dtype=np.int64)], [np.eye(n, dtype=np.int64), np.zeros((n, n), dtype=np.int64)]])
    return F12, FM, Om
def issymp(F, Om): return np.array_equal((F.T @ Om @ F) % 2, Om)
def inv2(A):
    n = A.shape[0]; M = np.concatenate([A % 2, np.eye(n, dtype=np.int64)], 1)
    r = 0
    for c in range(n):
        p = [i for i in range(r, n) if M[i, c]]
        if not p: raise ValueError
        M[[r, p[0]]] = M[[p[0], r]]
        for i in range(n):
            if i != r and M[i, c]: M[i] ^= M[r]
        r += 1
    return M[:, n:]
def ord3(A): I = np.eye(A.shape[0], dtype=np.int64); return not ((A @ A + A + I) % 2).any()
def rand_symp(dim, Om, rng, m=60):
    g = np.eye(dim, dtype=np.int64)
    for _ in range(m):
        v = rng.integers(0, 2, dim)
        T = (np.eye(dim, dtype=np.int64) + np.outer(v, v @ Om)) % 2   # x -> x + w(x,v) v ... as matrix: x + v (v^T Om x)
        g = (T @ g) % 2
    return g
if __name__ == '__main__':
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    F12, FM, Om = blocks(k); dim = 4*k
    print('symp', issymp(F12, Om), issymp(FM, Om), 'ord3 FM', ord3(FM), 'ord3 F12', ord3(F12))
    rng = np.random.default_rng(1); t0 = time.time(); hits = 0; tries = 0
    while time.time() - t0 < 120 and hits < 3:
        tries += 1
        g = rand_symp(dim, Om, rng); assert issymp(g, Om)
        A = (g @ FM @ inv2(g)) % 2
        B = (inv2(A) @ F12) % 2
        if ord3(B):
            hits += 1; print('hit', tries, 'A ord3', ord3(A), 'B symp', issymp(B, Om), 'A@B==F12', np.array_equal((A @ B) % 2, F12))
            np.save('A_k%d_%d.npy' % (k, hits), A); np.save('B_k%d_%d.npy' % (k, hits), B)
    print('tries', tries, 'hits', hits)
