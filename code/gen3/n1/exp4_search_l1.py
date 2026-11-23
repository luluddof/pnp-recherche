# Exhaustive search, lambda=1 (N=2): two-piece decompositions R + C = kappa*X, R in {-1,0,1}^16.
# Figure of merit: asymptotic rho of the blockwise tensor-power scheme (exact-rank eigen count),
# which reproduces 11/12 -> 0.901 (8-sq), 0.948 (4-sq), 0.896 (1-bit) of the paper.
import numpy as np, itertools, time, sys
from rho import *
N = 2
perm_row = np.array([[ (x*2+y)*4 + (v*2+u) for x in range(2) for y in range(2) for u in range(2) for v in range(2)]])
# index arithmetic: R flat index f = x*8+y*4+u*2+v ; build row/col matrices via fancy indexing
def mats(Rflat):
    n = Rflat.shape[0]
    Rt = Rflat.reshape(n, 2, 2, 2, 2)
    row = Rt.transpose(0, 1, 2, 4, 3).reshape(n, 4, 4)
    col = Rt.transpose(0, 1, 3, 4, 2).reshape(n, 4, 4)
    return row, col
vals = np.array([-1, 0, 1], dtype=np.float64)
seen = {}
t0 = time.time()
CH = 3**12
lowbits = np.array(list(itertools.product(range(3), repeat=12)), dtype=np.int8)
I4 = np.eye(4)
for kappa in (1, 2):
    for hi in itertools.product(range(3), repeat=4):
        digits = np.hstack([np.tile(np.array(hi, dtype=np.int8), (CH, 1)), lowbits])
        Rf = vals[digits]
        row, col = mats(Rf)
        rr = np.linalg.matrix_rank(row); rk = np.linalg.matrix_rank(row - kappa*I4)
        cr = np.linalg.matrix_rank(col); ck = np.linalg.matrix_rank(col - kappa*I4)
        key = np.stack([4-rr, 4-rk, 4-cr, 4-ck], 1)
        uk, first = np.unique(key, axis=0, return_index=True)
        for k_, f_ in zip(map(tuple, uk), first):
            if (kappa,)+k_ not in seen: seen[(kappa,)+k_] = Rf[f_].copy()
print(f"enumerated 2*3^16 candidates in {time.time()-t0:.0f}s; distinct spectral class tuples: {len(seen)}")
res = []
for key, Rf in seen.items():
    kappa = key[0]; R = Rf.reshape(2, 2, 2, 2)
    r, rc, cc, d, tau = rho_of(R, kappa, 1)
    res.append((r, kappa, rc, cc, d, tau, Rf))
res.sort(key=lambda z: z[0])
for r, kappa, rc, cc, d, tau, Rf in res[:8]:
    print(f"rho={r:.4f} kappa={kappa} row(z,k,m)={rc} col(z,k,m)={cc} diagonalizable={d} tau*={tau:.2f} R={Rf.astype(int).tolist()}")
np.save("l1_best.npy", res[0][6])
