# [E] exact checks for the "triangular swap lemma" (G24g notes sect. 2).
# Conventions: vector index (x,y) -> x*N+y, S = swap, K = sum_l kron(E_l, C_l), Theta = K S.
# Claim A: if all E_l are upper triangular, Theta restricted (principal) to Omega = {(b,j): b > j} is nilpotent,
#          so det(1 - Theta[Om,Om]) = 1 and rank(S - K) >= N(N-1)/2, for ANY C_l (any field).
# Claim B: every square pencil (E1,E2) in Kronecker canonical form becomes upper triangular after a row permutation
#          and a column permutation (order: L_eps blocks, regular blocks, L_eta^T blocks).
# Claim C (staircase): for any E_l, with t(b) = #leading zero columns common to row b of all E_l,
#          the principal minor on Omega' = {(b,j): j < min(b, t(b)+1)} (0-based: j < b, j < t(b)) is nilpotent.
import numpy as np, flint, itertools, random, sys

def S_mat(N):
    S = np.zeros((N*N, N*N), dtype=np.int64)
    for x in range(N):
        for y in range(N): S[y*N+x, x*N+y] = 1
    return S

def fm(A): return flint.fmpz_mat([[int(v) for v in row] for row in A])

def rank(A): return fm(A).rank()

def check_minor(Th, Om):
    idx = [b*N_+j for (b,j) in Om]
    sub = Th[np.ix_(idx, idx)]
    M = fm(sub)
    # nilpotent: M^n == 0
    P = M
    for _ in range(len(idx)): P = P*M
    nil = all(P[i,j] == 0 for i in range(len(idx)) for j in range(len(idx)))
    I = flint.fmpz_mat(len(idx), len(idx), [1 if i == j else 0 for i in range(len(idx)) for j in range(len(idx))])
    return nil, (I - M).det()

random.seed(1); rng = np.random.default_rng(1)
out = []
# ---- Claim A
for N_ in [2, 3, 4, 5]:
    N = N_; S = S_mat(N)
    for r in [1, 2, 3, 5]:
        for trial in range(6):
            Es = [np.triu(rng.integers(-2, 3, (N, N))) for _ in range(r)]
            Cs = [rng.integers(-2, 3, (N, N)) for _ in range(r)]
            K = sum(np.kron(E, C) for E, C in zip(Es, Cs))
            Th = K @ S
            Om = [(b, j) for b in range(N) for j in range(N) if b > j]
            nil, d = check_minor(Th, Om)
            rk = rank(S - K)
            ok = nil and d == 1 and rk >= N*(N-1)//2
            out.append(ok)
            if not ok: print("FAIL A", N, r, trial, nil, d, rk)
print("Claim A:", sum(out), "/", len(out), "ok")

# ---- also: adversarial C for upper-triangular E: try C = E^{-1}-like / identity choices
N_ = 4; N = 4; S = S_mat(N); best = 99
for trial in range(300):
    r = 3
    Es = [np.triu(rng.integers(-1, 2, (N, N))) for _ in range(r)]
    Cs = [rng.integers(-1, 2, (N, N)) for _ in range(r)]
    if trial % 3 == 0: Es[0] = np.eye(N, dtype=np.int64); Cs[0] = np.eye(N, dtype=np.int64)
    K = sum(np.kron(E, C) for E, C in zip(Es, Cs))
    best = min(best, rank(S - K))
print("Claim A random small-entry search, N=4, r=3, min rank found:", best, "(bound 6)")

# ---- Claim B: Kronecker blocks
def L(eps):   # eps x (eps+1)
    A = np.zeros((eps, eps+1), dtype=np.int64); B = np.zeros((eps, eps+1), dtype=np.int64)
    for i in range(eps): A[i, i] = 1; B[i, i+1] = 1
    return A, B
def LT(eta):  # (eta+1) x eta
    A, B = L(eta); return A.T.copy(), B.T.copy()
def J(k, lam):
    A = np.eye(k, dtype=np.int64); B = lam*np.eye(k, dtype=np.int64)
    for i in range(k-1): B[i, i+1] = 1
    return A, B
def Ninf(k):
    A, B = J(k, 0); return B, A
def blockdiag(blocks):
    R = sum(b[0].shape[0] for b in blocks); C = sum(b[0].shape[1] for b in blocks)
    A = np.zeros((R, C), dtype=np.int64); B = np.zeros((R, C), dtype=np.int64); r0 = c0 = 0
    for a, b in blocks:
        A[r0:r0+a.shape[0], c0:c0+a.shape[1]] = a; B[r0:r0+a.shape[0], c0:c0+a.shape[1]] = b
        r0 += a.shape[0]; c0 += a.shape[1]
    return A, B
okB = 0; totB = 0
for trial in range(400):
    nL = random.randint(0, 3)
    eps = [random.randint(0, 3) for _ in range(nL)]; eta = [random.randint(0, 3) for _ in range(nL)]
    reg = []
    for _ in range(random.randint(0, 3)):
        reg.append(J(random.randint(1, 3), random.randint(-2, 2)) if random.random() < .7 else Ninf(random.randint(1, 3)))
    blocks = [L(e) for e in eps] + reg + [LT(e) for e in eta]   # the claimed order
    random.shuffle(blocks) if False else None
    if not blocks: continue
    A, B = blockdiag(blocks)
    if A.shape[0] == 0: continue
    totB += 1
    assert A.shape[0] == A.shape[1]
    # block order as listed, natural row and column order inside: check upper triangular
    if np.all(np.tril(A, -1) == 0) and np.all(np.tril(B, -1) == 0): okB += 1
    else: print("FAIL B", eps, eta, [b[0].shape for b in reg])
print("Claim B:", okB, "/", totB, "pencils upper triangular in the stated block order")

# ---- Claim A on Kronecker pencils disguised by random equivalence (rank bound only)
okA2 = 0; tot = 0
for trial in range(40):
    eps = [random.randint(0, 2)]; eta = [random.randint(0, 2)]
    blocks = [L(eps[0])] + [J(random.randint(1, 2), random.randint(-1, 1))] + [LT(eta[0])]
    A, B = blockdiag(blocks); N = A.shape[0]
    if N > 5: continue
    N_ = N; S = S_mat(N)
    X = np.triu(rng.integers(-1, 2, (N, N)), 1) + np.eye(N, dtype=np.int64)
    Y = np.tril(rng.integers(-1, 2, (N, N)), -1) + np.eye(N, dtype=np.int64)
    E1, E2 = X @ A @ Y, X @ B @ Y
    C1, C2 = rng.integers(-2, 3, (N, N)), rng.integers(-2, 3, (N, N))
    K = np.kron(E1, C1) + np.kron(E2, C2)
    tot += 1; okA2 += rank(S - K) >= N*(N-1)//2
print("rank bound on disguised singular pencils:", okA2, "/", tot)

# ---- Claim C (staircase) on random general E_l
okC = 0; totC = 0
for N in [3, 4, 5]:
    N_ = N; S = S_mat(N)
    for trial in range(10):
        r = random.randint(2, 4)
        Es = [rng.integers(-1, 2, (N, N)) for _ in range(r)]
        # impose a staircase: row b of every E_l vanishes in columns < t(b)
        t = sorted(rng.integers(0, N, N))
        for E in Es:
            for b in range(N): E[b, :t[b]] = 0
        Cs = [rng.integers(-2, 3, (N, N)) for _ in range(r)]
        K = sum(np.kron(E, C) for E, C in zip(Es, Cs)); Th = K @ S
        tt = [min([c for c in range(N) if any(E[b, c] != 0 for E in Es)] + [N]) for b in range(N)]
        Om = [(b, j) for b in range(N) for j in range(N) if j < b and j < tt[b]]
        if not Om: continue
        nil, d = check_minor(Th, Om); rk = rank(S - K)
        totC += 1; okC += nil and d == 1 and rk >= len(Om)
print("Claim C:", okC, "/", totC)
