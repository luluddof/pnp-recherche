"""Grouped (polynomially many pieces) figure for decompositions at lambda0 = 1.

Two pieces A + B = kappa*T:  in direction i, X_i = M_i(A) M_i(T)^{-1} and kappa*I - X_i commute,
so the type-class piece G_j (j of k blocks use B) is the coefficient of t^j in
prod_blocks (X + t(kappa I - X)); if X_i is diagonalizable its rank is at most
#{eigen-tuples with possibly nonzero coefficient}  (exact unless accidental cancellation).
Exponent per block f_i(tau) = max { H(p) + sum p_c log w_c : p_z <= tau <= p_z + p_m }  (N1 method),
w = (z, k, m) = multiplicities of eigenvalue 0, kappa, other.
rho_grp = min_q max_tau sum_i q_i f_i(tau) / 4      (4 = trivial bits per block per direction).

General s pieces: multinomial upper bound  log rank_i(G_n) <= log multinom + sum_s n_s c_i(s),
rho_mult = min_q max_pi [ H(pi) + sum_s pi_s <q, c(s)> ] / 4  (H(pi) appears once per direction,
q sums to 1).  Plus exact ranks of G_n for k = 2 blocks as a sanity check."""
import itertools, math
import numpy as np
from scipy.optimize import linprog
from box3 import flat, rank, costs

G = np.linspace(0, 1, 201)
PZ, PM = np.meshgrid(G, G, indexing='ij'); OK = PZ + PM <= 1 + 1e-12
PK = np.clip(1 - PZ - PM, 0, 1)
def _H(p): return np.where(p > 0, -p * np.log2(np.where(p > 0, p, 1)), 0)
TAUS = np.linspace(0, 1, 101)

def f_curve(cls):
    z, k, m = cls
    val = _H(PZ) + _H(PM) + _H(PK)
    for p, w in ((PZ, z), (PK, k), (PM, m)):
        val = np.where(p > 0, np.where(w > 0, val + p * np.log2(max(w, 1)), -np.inf), val)
    val = np.where(OK, val, -np.inf)
    return np.array([np.max(np.where((PZ <= t + 1e-9) & (t <= PZ + PM + 1e-9), val, -np.inf))
                     for t in TAUS])

def classes(X, kappa, tol=1e-7):
    n = X.shape[0]
    z = n - rank(X); k = n - rank(X - kappa * np.eye(n))
    ev = np.linalg.eigvals(X)
    diag = int(np.sum(abs(ev) < 1e-6)) == z and int(np.sum(abs(ev - kappa) < 1e-6)) == k
    # full diagonalizability check
    w, V = np.linalg.eig(X)
    diag = diag and np.linalg.matrix_rank(V, tol=1e-6) == n
    return (z, k, n - z - k), diag

def minimax_rows(F):
    """F: (#candidates, 3) -> min_q max_row <q,F_row>."""
    F = np.asarray(F, float); n = F.shape[0]
    res = linprog([0, 0, 0, 1], A_ub=np.hstack([F, -np.ones((n, 1))]), b_ub=np.zeros(n),
                  A_eq=[[1, 1, 1, 0]], b_eq=[1], bounds=[(0, 1)] * 3 + [(None, None)])
    return res.fun, res.x[:3]

def rho_two_piece(A, T, kappa):
    fs, ok = [], True
    for i in range(3):
        P = flat(T, i)
        X = flat(A, i) @ np.linalg.inv(P)
        cls, d = classes(X, kappa)
        ok &= d
        fs.append(f_curve(cls))
    F = np.stack(fs, 1)                  # (taus, 3)
    v, q = minimax_rows(F)
    return v / 4, q, ok, F

def rho_multinomial(cvecs, steps=24):
    """min_q max_pi [H(pi) + sum pi_s <q,c_s>]/4 ; for fixed q the inner max is
    log2 sum_s 2^{<q,c_s>}  (Gibbs).  Minimise over a q-grid."""
    C = np.asarray(cvecs, float); best = (1e9, None)
    for a in range(steps + 1):
        for b in range(steps + 1 - a):
            q = np.array([a, b, steps - a - b]) / steps
            v = math.log2(np.sum(2.0 ** (C @ q)))
            if v < best[0]:
                best = (v, q)
    return best[0] / 4, best[1]
