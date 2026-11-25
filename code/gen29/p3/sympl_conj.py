# k=2: find A,B in Sp(8,2) with F12 = A B, A = g FM g^-1, and an explicit symplectic h with B = h FM h^-1.
import numpy as np, itertools, sys
sys.argv = ['x', '2', '5']
exec(open('sympl.py').read().split("assert inC(FM)")[0])
def nullspace_f2(Amat):
    A = Amat.copy() % 2; r, c = A.shape; piv = []; row = 0
    for col in range(c):
        p = next((i for i in range(row, r) if A[i, col]), None)
        if p is None: continue
        A[[row, p]] = A[[p, row]]
        for i in range(r):
            if i != row and A[i, col]: A[i] ^= A[row]
        piv.append(col); row += 1
        if row == r: break
    free = [j for j in range(c) if j not in piv]; basis = []
    for f in free:
        v = np.zeros(c, dtype=np.int64); v[f] = 1
        for i, pc_ in enumerate(piv):
            if A[i, f]: v[pc_] = 1
        basis.append(v)
    return basis
def conjugator(B):
    # solve h FM = B h  (h unknown n x n): linear in vec(h)
    rows = []
    for i in range(n):
        for j in range(n):
            # (h FM)[i,j] - (B h)[i,j] = sum_l h[i,l] FM[l,j] - sum_l B[i,l] h[l,j]
            v = np.zeros(n*n, dtype=np.int64)
            for l in range(n):
                v[i*n + l] ^= FM[l, j]
                v[l*n + j] ^= B[i, l]
            rows.append(v)
    basis = nullspace_f2(np.array(rows))
    print('solution space dim', len(basis))
    for t in range(200000):
        coeffs = rng.integers(0, 2, len(basis))
        h = (sum(c*b for c, b in zip(coeffs, basis)) % 2).reshape(n, n)
        if issymp(h): return h
    return None
for t in range(200000):
    g = rand_symp(); A = (g @ FM @ inv(g)) % 2; B = (inv(A) @ F12) % 2
    if inC(B): break
h = conjugator(B)
assert h is not None and issymp(h) and np.array_equal((h @ FM @ inv(h)) % 2, B)
assert np.array_equal((A @ B) % 2, F12)
print('VERIFIED k=2: F12 = (g FM g^-1)(h FM h^-1), g,h symplectic')
np.save('k2_g.npy', g); np.save('k2_h.npy', h)
