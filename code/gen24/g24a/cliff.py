# real antisymmetric anticommuting orthonormal Pauli matrices: max family of size rho(2^r)-1, by clique search
import numpy as np, random, functools
def rho(n):
    k = (n & -n).bit_length() - 1
    return [2*k+1, 2*k, 2*k, 2*k+2][k % 4]
def pauli_mat(ax, az, r):
    X = np.array([[0,1],[1,0]]); Z = np.array([[1,0],[0,-1]]); I = np.eye(2, dtype=int)
    M = np.array([[1]])
    for q in range(r):
        a = (ax >> q) & 1; b = (az >> q) & 1
        P = I
        if a: P = X
        if b: P = P @ Z
        M = np.kron(M, P)
    return M.astype(np.int64)
pc = lambda v: bin(v).count('1') & 1
@functools.lru_cache(None)
def family(r, seed=0):
    t = rho(2**r) - 1
    V = [(ax, az) for ax in range(2**r) for az in range(2**r) if pc(ax & az)]   # antisymmetric
    rng = random.Random(seed)
    def anti(a, b): return pc(a[0] & b[1]) ^ pc(a[1] & b[0])
    best = []
    for trial in range(2000):
        rng.shuffle(V)
        cl = []
        for v in V:
            if all(anti(v, w) for w in cl): cl.append(v)
        if len(cl) > len(best): best = cl
        if len(best) >= t: break
    assert len(best) >= t, (r, len(best), t)
    E = [pauli_mat(ax, az, r) for ax, az in best[:t]]
    n = 2**r
    for i, a in enumerate(E):
        assert (a.T == -a).all() and (a @ a.T == np.eye(n)).all()
        for b in E[:i]: assert (a @ b == -(b @ a)).all()
    return E
if __name__ == "__main__":
    for r in range(1, 9):
        E = family(r); print(r, 2**r, len(E), rho(2**r)-1)
