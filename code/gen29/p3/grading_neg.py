# Check L3' : for a ring C with banded node tensors (bandwidth w in grading f(t)=t), P_X C P_Y = 0 and
# rank(sigma - C) >= |Y|, Y = {t1 > max(t2,t3) + w}, X = sigma(Y).
import numpy as np, itertools
rng = np.random.default_rng(0)
N, w, b = 6, 1, 3
def banded():
    A = rng.standard_normal((N, N)); 
    for i in range(N):
        for j in range(N):
            if abs(i-j) > w+1: A[i, j] = 0
    return A
idx = lambda t: (t[0]*N + t[1])*N + t[2]
for trial in range(2):
    A = [[banded() for _ in range(b)] for _ in range(b)]  # A[l][i] : node1, bonds (31:l, 12:i)
    B = [[banded() for _ in range(b)] for _ in range(b)]  # B[i][j] : node2 (12:i, 23:j)
    D = [[banded() for _ in range(b)] for _ in range(b)]  # D[j][l] : node3 (23:j, 31:l)
    C = np.zeros((N**3, N**3))
    for l in range(b):
        for i in range(b):
            for j in range(b):
                C += np.kron(np.kron(A[l][i], B[i][j]), D[j][l])
    S = np.zeros((N**3, N**3))
    for t in itertools.product(range(N), repeat=3):
        S[idx((t[2], t[0], t[1])), idx(t)] = 1     # sigma|t1,t2,t3> = |t3,t1,t2>
    Y = [idx(t) for t in itertools.product(range(N), repeat=3) if t[0] > max(t[1], t[2]) + w]
    X = [idx((t[2], t[0], t[1])) for t in itertools.product(range(N), repeat=3) if t[0] > max(t[1], t[2]) + w]
    leak = np.abs(C[np.ix_(X, Y)]).max()
    for scale in (1.0, 0.3):
        r = np.linalg.matrix_rank(S - scale*C)
        print(f'trial {trial}: |Y|={len(Y)} leak P_X C P_Y max={leak:.1e} rank(sigma - {scale}C)={r} (>= |Y| required)')
# negative control: bandwidth w+1 factors should leak
