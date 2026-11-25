# achievable functionals A_k (coefficient vectors in the basis e_nu^*) and the "bad" directions b in V_h with l(b)=0 for all l in A_k
import sys, pickle; sys.dont_write_bytecode = True
import numpy as np, flint
from feas import allowed
from appb_exact import partitions
fn, W, Wp = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
d = pickle.load(open(fn, 'rb')); h = d['h']; P = d['P']; Vh = sorted(partitions(h)); n = len(Vh)
cols = [v for (mu, e), v in d['cols'].items() if allowed(h, mu, e, W, Wp)]
A = np.stack(cols + d['rhs'], axis=1) % P
M = flint.nmod_mat(A.shape[0], A.shape[1], A.ravel().tolist(), P)
X, nul = M.nullspace()
Y = np.array([[int(X[len(cols) + i, j]) for i in range(n)] for j in range(nul)], dtype=np.int64)   # each row: an achievable l
My = flint.nmod_mat(Y.shape[0], n, (Y % P).ravel().tolist(), P)
r = My.rank(); print(f"h={h} W={W} W'={Wp} k={4*h-1-W-Wp}: dim A_k = {r} of {n}")
Z, nz = My.nullspace()      # b with Y b = 0
for j in range(nz):
    b = [int(Z[i, j]) for i in range(n)]
    b = [x if x <= P // 2 else x - P for x in b]
    print("  bad b:", {Vh[i]: b[i] for i in range(n) if b[i]})
