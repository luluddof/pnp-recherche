# Kronecker square of the exact (5,3) identity at N = 4 (N = 16): exact ranks of the pieces A_i (i = number of C factors) and thresholds.
import sys, numpy as np
sys.path.insert(0, '/home/user/pnp-recherche/code/gen18/g18g')
import exact
G = [[[0,0,0,0],[0,0,0,0],[0,0,-1,0],[0,1,0,0]], [[0,0,0,-1],[0,0,-1,0],[0,0,0,0],[1,0,0,0]], [[0,0,0,0],[0,1,0,0],[-1,0,0,0],[0,0,0,0]]]
Gp = [[[0,0,0,0],[0,0,0,1],[0,0,-1,0],[0,0,0,0]], [[0,0,0,1],[0,0,0,0],[0,-1,0,0],[-1,0,0,0]], [[0,0,-1,0],[0,1,0,0],[0,0,0,0],[0,0,0,0]]]
C = np.einsum('kxu,kyv->xyuv', np.array(G), np.array(Gp)).astype(np.int64); R = exact.X(4) - C
print("base: rr(R), rc(R) =", exact.rank(exact.rowm(R)), exact.rank(exact.colm(R)), " rr(C), rc(C) =", exact.rank(exact.rowm(C)), exact.rank(exact.colm(C)))
exact.power("(5,3)@N=4", (R, C, 1), 1); exact.power("(5,3)@N=4", (R, C, 1), 2)
