"""lambda = 1: SU(2) skein view.  Local invertible maps on single corners preserve every
flattening rank, so replace T by T' = prod_diag eps_{p,7-p} (apply eps to the R_1 endpoint).
SL2-invariants on 8 qubits = span of eps-matchings (dim Catalan(4)=14).  An eps-matching
eps_M has flattening rank 2^{cross_i(M)} in direction i.  Find representations
T' = sum_M a_M eps_M and their single-block ('exponential expansion') figure
    rho_exp = min_q max_{M in supp} <q, cross(M)> / 4."""
import itertools, math
import numpy as np
from scipy.optimize import linprog
from box3 import bit, all_matchings, costs, flat, T_tensor, CORNERS

EPS = np.array([[0., 1.], [-1., 0.]])

def pair_tensor(p, q, Mx):
    sh = [1] * 8; sh[p] = 2; sh[q] = 2
    return (Mx if p < q else Mx.T).reshape(sh)

def eps_matching(M):
    E = np.ones((2,) * 8)
    for p, q in M:
        E = E * pair_tensor(min(p, q), max(p, q), EPS)
    return E

def cross(M):
    return tuple(sum(bit(p, i) != bit(q, i) for p, q in M) for i in range(3))

DIAG = [(p, 7 - p) for p in range(4)]
Tp = eps_matching(DIAG)              # = T with eps applied on corners 4..7 (up to sign)

def in_span(target, mats):
    if not mats:
        return False, None
    A = np.stack([m.ravel() for m in mats], 1)
    x, *_ = np.linalg.lstsq(A, target.ravel(), rcond=None)
    return np.allclose(A @ x, target.ravel(), atol=1e-9), x

def minimax(vecs):
    """min_q max_v <q,v>, q in simplex (LP)."""
    V = np.array(vecs, float); n = V.shape[0]
    # vars q1,q2,q3,z ; minimize z ; V q - z <= 0 ; sum q = 1
    res = linprog([0, 0, 0, 1], A_ub=np.hstack([V, -np.ones((n, 1))]), b_ub=np.zeros(n),
                  A_eq=[[1, 1, 1, 0]], b_eq=[1], bounds=[(0, 1)] * 3 + [(None, None)])
    return res.fun, res.x[:3]

if __name__ == "__main__":
    Ms = list(all_matchings(CORNERS))
    print("check: costs(T') =", costs(Tp), " costs(T) =", costs(T_tensor(1)))
    for M in Ms[:3]:
        assert [round(c) for c in costs(eps_matching(M))] == list(cross(M))
    # dimension of span of all eps-matchings
    A = np.stack([eps_matching(M).ravel() for M in Ms], 1)
    print("dim span eps-matchings =", np.linalg.matrix_rank(A))
    # threshold on max_i cross_i and on sum
    for label, key in [("sum", lambda c: sum(c)), ("max", lambda c: max(c))]:
        for thr in range(0, 13):
            S = [M for M in Ms if key(cross(M)) <= thr]
            ok, x = in_span(Tp, [eps_matching(M) for M in S])
            if ok:
                print(f"T' in span of matchings with {label}(cross) <= {thr}: |S|={len(S)}")
                break
    # cube-edge matchings
    edges = [M for M in Ms if sum(cross(M)) == 4]
    print("#cube-edge perfect matchings:", len(edges), sorted(set(cross(M) for M in edges)))
    ok, x = in_span(Tp, [eps_matching(M) for M in edges])
    print("T' in span of cube-edge matchings?", ok)
    # best rho_exp over supports: greedy search over subsets defined by q-weighted threshold
    best = (9, None)
    qs = [np.array(q) / sum(q) for q in itertools.product(range(0, 7), repeat=3) if sum(q) > 0]
    for q in qs:
        vals = sorted(set(round(float(np.dot(q, cross(M))), 9) for M in Ms))
        for v in vals:
            S = [M for M in Ms if np.dot(q, cross(M)) <= v + 1e-9]
            ok, x = in_span(Tp, [eps_matching(M) for M in S])
            if ok:
                # actual support used: minimal by LP over this set (use lstsq support)
                r, qq = minimax([cross(M) for M in S])
                if r < best[0] - 1e-9:
                    best = (r, q, len(S), sorted(set(cross(M) for M in S)))
                break
    print("best single-block figure min_q max <q,c>/4 =", best[0] / 4, "q=", best[1],
          "|S|=", best[2], "cost vectors", best[3])
