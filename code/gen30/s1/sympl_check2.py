import numpy as np
from sympl_check import *
def nullspace2(M):
    M = M.copy() % 2; r = 0; piv = []; m, n = M.shape
    for c in range(n):
        p = [i for i in range(r, m) if M[i, c]]
        if not p: continue
        M[[r, p[0]]] = M[[p[0], r]]
        for i in range(m):
            if i != r and M[i, c]: M[i] ^= M[r]
        piv.append(c); r += 1
    free = [c for c in range(n) if c not in piv]; basis = []
    for f in free:
        v = np.zeros(n, dtype=np.int64); v[f] = 1
        for i, c in enumerate(piv): v[c] = M[i, f]
        basis.append(v)
    return basis
def conj_to(FM, B, Om, rng, tries=20000):
    d = FM.shape[0]
    # unknown X (d x d), equations X FM - B X = 0, vec row-major
    rows = []
    for i in range(d):
        for j in range(d):
            r = np.zeros(d*d, dtype=np.int64)
            for l in range(d): r[i*d + l] ^= FM[l, j]      # (X FM)_{ij} = sum_l X_il FM_lj
            for l in range(d): r[l*d + j] ^= B[i, l]       # (B X)_{ij} = sum_l B_il X_lj
            rows.append(r)
    basis = nullspace2(np.array(rows))
    for _ in range(tries):
        c = rng.integers(0, 2, len(basis)); X = (c @ np.array(basis)) % 2; X = X.reshape(d, d)
        if issymp(X, Om): return X
    return None
rng = np.random.default_rng(5)
F12, FM, Om = blocks(2)
for h in (1, 2, 3):
    A = np.load('A_k2_%d.npy' % h); B = np.load('B_k2_%d.npy' % h)
    X = conj_to(FM, B, Om, rng); Y = conj_to(FM, A, Om, rng)
    print(h, 'B = X FM X^-1 with X symplectic:', X is not None and np.array_equal((X @ FM @ inv2(X)) % 2, B),
          '| A likewise:', Y is not None and np.array_equal((Y @ FM @ inv2(Y)) % 2, A))
# k=4 by pairing sites: k=4 coords (ax0..3, ay0..3, bx0..3, by0..3); sites {0,1} and {2,3}
A = np.load('A_k2_1.npy'); B = np.load('B_k2_1.npy')
F12_4, FM_4, Om4 = blocks(4)
def embed(A2, sites):
    d4 = 16; E = np.eye(d4, dtype=np.int64)
    idx4 = [blk*4 + s for blk in range(4) for s in sites]    # k=2 coordinate (blk, j) -> k=4 (blk, sites[j])
    E[np.ix_(idx4, idx4)] = A2; return E
A4 = (embed(A, [0, 1]) @ embed(A, [2, 3])) % 2; B4 = (embed(B, [0, 1]) @ embed(B, [2, 3])) % 2
print('k=4: A4 B4 = F12', np.array_equal((A4 @ B4) % 2, F12_4), 'ord3', ord3(A4), ord3(B4), 'symp', issymp(A4, Om4), issymp(B4, Om4))
X = conj_to(FM_4, A4, Om4, rng, 50000); print('k=4 A4 conj to FM (explicit symplectic X):', X is not None)
