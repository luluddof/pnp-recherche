# [E] Exhaustive check of Theorem 1 (r = 2) at N = 3 over GF(2) (characteristic 2: G19b's r = 1 proof does not apply).
# Every pencil (E1,E2) is equivalent (X E Y^{-1}) to one with E1 = diag(I_k,0); E2 arbitrary (512); C1, C2 arbitrary (2^18).
# Computes min over all of rank_{GF(2)}(Swap - E1(x)C1 - E2(x)C2).  Theorem 1 predicts exactly 3.
import numpy as np, itertools, time
N = 3; M = N * N
def swap():
    S = np.zeros((M, M), dtype=np.uint8)
    for x in range(N):
        for y in range(N): S[x*N+y, y*N+x] = 1
    return S
def rows_to_int(A): return np.array([int(''.join(map(str, r)), 2) for r in A], dtype=np.uint16)
def unit(i, j):
    U = np.zeros((N, N), dtype=np.uint8); U[i, j] = 1; return U
units = [unit(i, j) for i in range(N) for j in range(N)]
def batch_rank(R):   # R: (B, M) uint16 rows as bitmasks over M columns
    R = R.copy(); B = R.shape[0]; rank = np.zeros(B, dtype=np.int32); used = np.zeros((B, M), dtype=bool)
    for bit in range(M):
        mask = np.uint16(1 << bit)
        has = ((R & mask) != 0) & ~used
        any_ = has.any(axis=1)
        piv = np.argmax(has, axis=1)
        idx = np.arange(B)
        prow = R[idx, piv]
        elim = ((R & mask) != 0) & any_[:, None]
        elim[idx, piv] = False
        R = np.where(elim, R ^ prow[:, None], R)
        used[idx[any_], piv[any_]] = True
        rank += any_
    return rank
S = swap(); Sint = rows_to_int(S)
best = 99; t0 = time.time(); count = 0
c = np.arange(1 << 18, dtype=np.int64)
allm = [np.array([(v >> i) & 1 for i in range(9)], dtype=np.uint8).reshape(N, N) for v in range(512)]
def code(A): return int(sum(int(A.ravel()[i]) << i for i in range(9)))
def inv2(A):
    import flint
    m = flint.nmod_mat([[int(v) for v in r] for r in A], 2)
    if m.rank() < N: return None
    mi = m.inv(); return np.array([[int(mi[i, j]) for j in range(N)] for i in range(N)], dtype=np.uint8)
GL = [(A, inv2(A)) for A in allm if inv2(A) is not None]
print("|GL3(F2)| =", len(GL), flush=True)
reps = []
for k in range(N + 1):
    D = np.diag([1]*k + [0]*(N-k)).astype(np.uint8)
    stab = [(X, Yi) for X, Xi in GL for Y, Yi in GL if np.array_equal((X @ D) % 2, (D @ Y) % 2)]
    seen = set()
    for v in range(512):
        if v in seen: continue
        reps.append((k, v)); orbit = {v}; frontier = [v]
        while frontier:
            w = frontier.pop(); W = allm[w]
            nb = [code((X @ W @ Yi) % 2) for X, Yi in stab] + [code((W + D) % 2)]
            for u in nb:
                if u not in orbit: orbit.add(u); frontier.append(u)
        seen |= orbit
print("pencil orbit representatives:", len(reps), flush=True)
for k, e2 in reps:
        E1 = np.diag([1]*k + [0]*(N-k)).astype(np.uint8)
        E2 = allm[e2]
        basis = [np.kron(E1, U) % 2 for U in units] + [np.kron(E2, U) % 2 for U in units]   # 18 generators
        bint = np.array([rows_to_int(Bm) for Bm in basis])   # (18, M)
        # rows = S ^ XOR_{i: c_i = 1} basis_i
        T1 = np.zeros((512, M), dtype=np.uint16); T2 = np.zeros((512, M), dtype=np.uint16)
        for s in range(512):
            for i in range(9):
                if (s >> i) & 1: T1[s] ^= bint[i]; T2[s] ^= bint[9+i]
        R = T1[c & 511] ^ T2[c >> 9] ^ Sint[None, :]
        rk = batch_rank(R).min(); count += 1
        if rk < best: best = rk; print(f"k={k} e2={e2} new min rank {rk}", flush=True)
print(f"N=3 GF(2) r=2: min rank over all {count} normalised pencils x 2^18 C-pairs = {best}  ({time.time()-t0:.0f}s)")
