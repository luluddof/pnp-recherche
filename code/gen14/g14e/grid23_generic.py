#!/usr/bin/env python3
"""G14e-3: 2x3 grid, lam=1, generic (non-matching) search.  T = A + B, rk_row(A) <= r,
B = ONE triangle network over the 3 columns with bond dims (D12,D13,D23).  Alternating: A <- truncated SVD of
rowflat(T-B); B <- ALS sweeps on T-A.  Also: grouped ranks of the product (one-bit identity) expansion."""
import sys, itertools, numpy as np
rng = np.random.default_rng(0)
# axes: t1,b1,t2,b2,t3,b3 each of dim 4 = 2 bits; node (row a, col j) bit k <-> partner column = the k-th other column
T = np.zeros((4,)*6)
oth = {0: [1, 2], 1: [0, 2], 2: [0, 1]}
for bits in itertools.product([0, 1], repeat=6):          # one bit per diagonal edge (top j -> bottom j'), j != j'
    e = dict(zip([(j, j2) for j in range(3) for j2 in oth[j]], bits))
    t = [sum(e[(j, j2)] << k for k, j2 in enumerate(oth[j])) for j in range(3)]
    b = [sum(e[(j2, j)] << k for k, j2 in enumerate(oth[j])) for j in range(3)]
    T[t[0], b[0], t[1], b[1], t[2], b[2]] = 1
rowflat = lambda X: X.transpose(0, 2, 4, 1, 3, 5).reshape(64, 64)
unrow = lambda M: M.reshape(4, 4, 4, 4, 4, 4).transpose(0, 3, 1, 4, 2, 5)
col3 = lambda X: X.reshape(16, 16, 16)
def cutranks(X):
    X = col3(X); return [np.linalg.matrix_rank(np.moveaxis(X, j, 0).reshape(16, 256), tol=1e-8) for j in range(3)]
print("T: row rank", np.linalg.matrix_rank(rowflat(T)), "col cut ranks", cutranks(T))
def tn(f1, f2, f3): return np.einsum('aij,bik,cjk->abc', f1, f2, f3)
def search(r, D, tries=6, iters=400):
    D12, D13, D23 = D; best = 1e9
    for t in range(tries):
        f1 = rng.normal(size=(16, D12, D13)); f2 = rng.normal(size=(16, D12, D23)); f3 = rng.normal(size=(16, D13, D23))
        A = np.zeros((4,)*6)
        for it in range(iters):
            Rm = col3(T - A)
            for _ in range(2):
                E = np.einsum('bik,cjk->ijbc', f2, f3).reshape(D12*D13, 256)
                f1 = np.linalg.lstsq(E.T, Rm.reshape(16, 256).T, rcond=None)[0].T.reshape(16, D12, D13)
                E = np.einsum('aij,cjk->ikac', f1, f3).reshape(D12*D23, 256)
                f2 = np.linalg.lstsq(E.T, Rm.transpose(1, 0, 2).reshape(16, 256).T, rcond=None)[0].T.reshape(16, D12, D23)
                E = np.einsum('aij,bik->jkab', f1, f2).reshape(D13*D23, 256)
                f3 = np.linalg.lstsq(E.T, Rm.transpose(2, 0, 1).reshape(16, 256).T, rcond=None)[0].T.reshape(16, D13, D23)
            B = tn(f1, f2, f3).reshape((4,)*6)
            if r > 0:
                u, s, vt = np.linalg.svd(rowflat(T - B)); A = unrow((u[:, :r]*s[:r]) @ vt[:r])
            err = np.sum((T - A - B)**2)
            if err < 1e-12: break
        best = min(best, err)
        if best < 1e-12: break
    return best
if __name__ == "__main__":
    # product baseline, grouped:  T = prod_gamma (a_gamma - b_gamma)  in delta/eps form, checked numerically
    cases = eval(sys.argv[1]) if len(sys.argv) > 1 else [(64, (1, 1, 1)), (0, (4, 4, 4))]
    for r, D in cases:
        e = search(r, D)
        print(f"r={r:2d} bonds={D} cost r+prod(D)={r+D[0]*D[1]*D[2]:3d}: best residual^2={e:.3e}", "FEASIBLE" if e < 1e-10 else "not found", flush=True)
