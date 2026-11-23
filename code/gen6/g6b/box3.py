"""G6b: 3-D crossing ('box') T on the 2x2x2 cube, lambda = 1 (and general helpers).
Corner p in {0,1}^3 is encoded as integer p = 4*p1 + 2*p2 + p3.  A function A on
({0,1}^lam)^8 is stored as an array of shape (N,)*8, N = 2^lam, axis p = corner p.
T(w) = prod_{p: p1=0} [w_p = w_{7-p}]   (4 main diagonals {p, 7-p}).
Direction i (merge along axis i): bipartition L_i = {p: p_i = 0} | R_i = {p: p_i = 1}.
c_i(A) = log2 rank_Q of the flattening M_i(A) (rows: corners of L_i, cols: corners of R_i,
columns ordered by the complement of the row corner, so that M_i(T) = identity)."""
import itertools, math
import numpy as np

CORNERS = list(range(8))
def bit(p, i):  # i in {0,1,2} -> coordinate i+1
    return (p >> (2 - i)) & 1
L = [[p for p in CORNERS if bit(p, i) == 0] for i in range(3)]
# column order: complement of each row corner, so T flattens to the identity
R = [[7 - p for p in L[i]] for i in range(3)]

def T_tensor(lam=1):
    N = 2 ** lam
    T = np.zeros((N,) * 8)
    for vals in itertools.product(range(N), repeat=4):
        idx = [0] * 8
        for p, v in zip([0, 1, 2, 3], vals):
            idx[p] = v; idx[7 - p] = v
        T[tuple(idx)] = 1
    return T

def flat(A, i):
    N = A.shape[0]
    return A.transpose(L[i] + R[i]).reshape(N ** 4, N ** 4)

def rank(M, tol=1e-7):
    return int(np.linalg.matrix_rank(M, tol=tol))

def costs(A):
    return [math.log2(max(rank(flat(A, i)), 1)) for i in range(3)]

def matching_fn(M, lam=1):
    """E_M = prod_{(p,q) in M} [w_p = w_q]."""
    N = 2 ** lam
    E = np.ones((N,) * 8)
    for p, q in M:
        sh = [1] * 8
        eq = np.eye(N)
        sh[p] = N; sh[q] = N
        # build [w_p = w_q] broadcast
        a = np.zeros((N, N)); a[:] = eq
        if p < q:
            E = E * a.reshape([N if k in (p, q) else 1 for k in range(8)])
        else:
            E = E * a.T.reshape([N if k in (p, q) else 1 for k in range(8)])
    return E

def all_matchings(pts):
    if not pts:
        yield []
        return
    a = pts[0]
    for j in range(1, len(pts)):
        b = pts[j]
        rest = pts[1:j] + pts[j + 1:]
        for m in all_matchings(rest):
            yield [(a, b)] + m

if __name__ == "__main__":
    T = T_tensor(1)
    print("T costs", costs(T), "flat is identity:", all(np.allclose(flat(T, i), np.eye(16)) for i in range(3)))
    Ms = list(all_matchings(CORNERS))
    print("#matchings", len(Ms))
    from collections import Counter
    cnt = Counter(tuple(int(round(c)) for c in costs(matching_fn(M))) for M in Ms)
    print("cost vectors of matching products E_M:", sorted(cnt.items()))
    # 2-D analogue check: T - E_axis1
    E1 = matching_fn([(p, p + 4) for p in range(4)])
    print("E_axis1 costs", costs(E1), " T-E_axis1 costs", costs(T - E1))
