"""Exact check: tensor square of the one-bit identity X = R - C', R = [x=u][y=v], C' = g(x,y)g(u,v), g = [x!=y](-1)^x.
N = 4: the 4 ungrouped pieces and the grouped middle piece; ranks over Q (sympy) and GF(2)."""
import itertools, numpy as np, sympy
from sat_k import gf2rank
R = np.zeros((2,)*4, dtype=np.int64); C = np.zeros((2,)*4, dtype=np.int64); X1 = np.zeros((2,)*4, dtype=np.int64)
g = lambda a, b: (a != b) * (-1) ** a
for x, y, u, v in itertools.product(range(2), repeat=4):
    R[x, y, u, v] = (x == u) * (y == v); C[x, y, u, v] = -g(x, y) * g(u, v); X1[x, y, u, v] = (x == v) * (y == u)
assert (R + C == X1).all()
def kron(A, B):   # indices (x1x2, y1y2, u1u2, v1v2)
    T = np.einsum('abcd,efgh->aebfcgdh', A, B).reshape(4, 4, 4, 4); return T
def ranks(T):
    row = T.reshape(16, 16); col = T.transpose(0, 2, 1, 3).reshape(16, 16)
    return (sympy.Matrix(row.tolist()).rank(), sympy.Matrix(col.tolist()).rank(), gf2rank(row % 2), gf2rank(col % 2))
X = kron(X1, X1); ps = {'RR': kron(R, R), 'RC': kron(R, C), 'CR': kron(C, R), 'CC': kron(C, C)}
assert (sum(ps.values()) == X).all()
Xc = np.zeros((4,)*4, dtype=np.int64)
for x, y in itertools.product(range(4), repeat=2): Xc[x, y, y, x] = 1
assert (X == Xc).all()
for k, T in ps.items(): print(k, "ranks (row_Q, col_Q, row_GF2, col_GF2) =", ranks(T))
print("grouped middle RC+CR:", ranks(ps['RC'] + ps['CR']))
