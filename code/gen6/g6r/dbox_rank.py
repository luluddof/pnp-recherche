# c_i(T) for the d-box crossing tensor T(w)=prod_{diagonals {p,~p}} [w_p = w_~p], corners in {0,1}^d, lambda bits each
import itertools, numpy as np
def T_tensor(d, lam):
    N = 2**lam; C = 2**d; full = C-1
    T = np.zeros((N,)*C)
    reps = [p for p in range(C) if p < full-p]
    for vals in itertools.product(range(N), repeat=len(reps)):
        idx = [0]*C
        for p, v in zip(reps, vals): idx[p] = v; idx[full-p] = v
        T[tuple(idx)] = 1
    return T
def c(T, d, i):
    C = 2**d; N = T.shape[0]
    L = [p for p in range(C) if not (p >> i) & 1]; R = [p for p in range(C) if (p >> i) & 1]
    M = T.transpose(L+R).reshape(N**len(L), N**len(R))
    return np.log2(np.linalg.matrix_rank(M))
for d, lam in [(2,1),(2,2),(3,1),(3,2),(4,1)]:
    T = T_tensor(d, lam)
    print(f"d={d} lam={lam}: c_i(T) =", [round(c(T,d,i),3) for i in range(d)], " 2^(d-1)*lam =", 2**(d-1)*lam, " d*lam =", d*lam)
# minimax sanity: random nonneg cost vectors where each direction i has some piece with c_i>=L
from scipy.optimize import linprog
rng = np.random.default_rng(0); worst = 9
for _ in range(2000):
    d = rng.integers(2,6); L = 1.0
    V = rng.random((rng.integers(d, 3*d), d)) * rng.random()  # small random costs
    for i in range(d): V[i, :] = rng.random(d)*0.1; V[i, i] = L   # piece i large only in direction i
    k = V.shape[0]
    res = linprog(np.r_[np.zeros(d),1], A_ub=np.c_[V, -np.ones(k)], b_ub=np.zeros(k),
                  A_eq=np.r_[np.ones(d),0][None], b_eq=[1], bounds=[(0,None)]*d+[(None,None)])
    worst = min(worst, res.fun * d / L)
print("min over trials of d*minimax/L =", round(worst,4), "(claim >= 1)")
