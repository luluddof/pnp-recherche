# Load / build pieces SWAP' = R + C  (R: low rank in row flattening (x,y)|(u,v); C: low osr = rank in column flattening (x,u)|(y,v)).
# Tensors indexed [x,y,u,v]; Swap[x,y,u,v] = [x=v][y=u].
import json, numpy as np, itertools
def swap(N):
    S = np.zeros((N,)*4, dtype=np.int64)
    for x in range(N):
        for y in range(N): S[x, y, y, x] = 1
    return S
def from_EC(N, E, C, scale=1):
    K = np.einsum('kxu,kyv->xyuv', np.array(E, dtype=np.int64), np.array(C, dtype=np.int64))
    return scale*swap(N) - K, K, scale
def load_json(fn):
    d = json.load(open(fn)); return d['N'], from_EC(d['N'], d['E'], d['C'])
def F1(N):   # (N(N-1)/2, 1): K = identity operator
    I = np.eye(N, dtype=np.int64); return N, from_EC(N, [I], [I])
def octonion():   # Lemma 2 / (30): R + C = 2 Swap on blocks of 3 bits (N = 8)
    def phi(p, q):
        s = sum(p[i]*q[j] for i in range(3) for j in range(i, 3)) + q[0]*p[1]*p[2] + p[0]*q[1]*p[2] + p[0]*p[1]*q[2]
        return (-1)**(s % 2)
    B = list(itertools.product([0, 1], repeat=3)); idx = {b: i for i, b in enumerate(B)}
    xor = lambda a, b: tuple(i ^ j for i, j in zip(a, b))
    R = np.zeros((8,)*4, dtype=np.int64); C = np.zeros((8,)*4, dtype=np.int64)
    for x, y, u, v in itertools.product(B, repeat=4):
        i = (idx[x], idx[y], idx[u], idx[v])
        if xor(x, y) == xor(v, u): R[i] = phi(x, y)*phi(v, u)
        if xor(x, u) == xor(v, y): C[i] = phi(x, u)*phi(v, y)
    assert (R + C == 2*swap(8)).all()
    return 8, (R, C, 2)
def rowflat(T): N = T.shape[0]; return T.reshape(N*N, N*N)                       # (x,y)|(u,v)
def colflat(T): N = T.shape[0]; return T.transpose(0, 2, 1, 3).reshape(N*N, N*N)  # (x,u)|(y,v)
