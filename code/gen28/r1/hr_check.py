# Hurwitz-Radon tiling: integer SOS [N,N,s], s = N*ceil(N/rho(N)), and the induced crossing pair (s,s):
# D = sum_k vec(F_k) vec(F_k)^T (rank <= s); R = D*SWAP/2 ; C = SWAP - R ; check rank R, osr C exactly.
import numpy as np, itertools, sys, time
from math import ceil
from singlet_check import rank_mod, realign, PR
from kraw import rho
X = np.array([[0,1],[1,0]]); Z = np.array([[1,0],[0,-1]]); Y = X @ Z; I2 = np.eye(2, dtype=int)
P = {'I': I2, 'X': X, 'Z': Z, 'Y': Y}
def mat(s):
    M = np.array([[1]])
    for ch in s: M = np.kron(M, P[ch])
    return M
def anti(s, t): return sum(1 for a, b in zip(s, t) if a != 'I' and b != 'I' and a != b) % 2 == 1
def hr_family(k, need):
    sk = [''.join(t) for t in itertools.product('IXZY', repeat=k) if ''.join(t).count('Y') % 2 == 1]
    rng = np.random.default_rng(1)
    for trial in range(2000):   # randomized greedy with restarts
        order = rng.permutation(len(sk)); chosen = []
        for i in order:
            s = sk[i]
            if all(anti(s, t) for t in chosen): chosen.append(s)
            if len(chosen) == need: return chosen
    return None
def build(k):
    N = 2**k; r = rho(N)
    gens = hr_family(k, r - 1); assert gens is not None, "no HR family found"
    A = [np.eye(N, dtype=np.int64)] + [mat(s).astype(np.int64) for s in gens]
    for i in range(r):
        for j in range(r):
            S = A[i].T @ A[j] + A[j].T @ A[i]
            assert (S == (2*np.eye(N, dtype=np.int64) if i == j else 0)).all()
    # forms: block b of rows (x indices) of size <= r ; z_{b,l} = sum_{i in block} x_i (A_{i-b0} y)_l
    F = []
    for b0 in range(0, N, r):
        for l in range(N):
            Fm = np.zeros((N, N), dtype=np.int64)
            for i in range(b0, min(b0 + r, N)): Fm[i, :] = A[i - b0][l, :]
            F.append(Fm.reshape(-1))
    Lm = np.array(F)   # s x N^2
    return N, r, Lm
if __name__ == "__main__":
    k = int(sys.argv[1]); full = len(sys.argv) > 2
    t = time.time(); N, r, Lm = build(k); s = Lm.shape[0]
    print(f"N={N} rho={r} HR conditions OK; integer SOS [N,N,{s}] (N^2/8={N*N//8})", flush=True)
    # SOS identity check: sum_k (x^T F_k y)^2 == |x|^2|y|^2 on random integer vectors
    rng = np.random.default_rng(0)
    for _ in range(3):
        x = rng.integers(-5, 6, N); y = rng.integers(-5, 6, N)
        v = Lm @ np.kron(x, y); assert (v @ v) == (x @ x) * (y @ y)
    print("SOS identity holds on random vectors", flush=True)
    if full:
        n = N*N
        Lf = Lm.astype(np.float64); D = np.rint(Lf.T @ Lf).astype(np.int64)
        x = np.repeat(np.arange(N), N); y = np.tile(np.arange(N), N)
        SW = np.zeros((n, n), dtype=np.int64); SW[y*N + x, np.arange(n)] = 1
        R2 = D @ SW if N <= 16 else D[:, (y*N + x)]   # 2R = D*SWAP (column permutation)
        C2 = 2*SW - R2       # 2C
        for p in PR:
            print(f"mod {p}: rank R = {rank_mod(R2 % p, p)}, osr C = {rank_mod(realign(C2, N) % p, p)}", flush=True)
    print("time %.1fs" % (time.time() - t))
