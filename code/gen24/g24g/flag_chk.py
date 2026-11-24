# [E] Lemma C (flag form): for E-span with structured flags, build X from adapted bases of (E F_t)^perp, check that the
# principal minor of Theta = K Swap on Om' is nilpotent and rank(Swap - K) >= sum_t (N - max(t, dim E F_t)).
# Families tested: E_l = A O_l B with O_l upper triangular on a module of dim m > N (A: m->N, B: N->m random), so that
# dim E F_t < r t is possible; C_l random integers.
import numpy as np, flint, sympy as sp, random
rng = np.random.default_rng(3)
def Smat(N):
    S = np.zeros((N*N, N*N), dtype=np.int64)
    for x in range(N):
        for y in range(N): S[y*N+x, x*N+y] = 1
    return S
def fm(A): return flint.fmpq_mat([[sp.Rational(v).p if False else int(v) for v in row] for row in A]) if False else None
def qrank(A): return sp.Matrix(A).rank()
ok = tot = 0
for trial in range(8):
    N = random.choice([3, 4]); m = N + random.choice([0, 1, 2]); r = random.choice([2, 3, 4])
    O = [np.triu(rng.integers(-2, 3, (m, m))) for _ in range(r)]
    A = rng.integers(-2, 3, (N, m)); B = rng.integers(-2, 3, (m, N))
    Es = [A @ o @ B for o in O]
    # flag F_t = span of first t columns of Y^{-1}; take Y^{-1} = identity after choosing B-adapted basis: use F_t = B^+-preimage
    # simple choice: F_t = span(e_1..e_t) in a random basis Fb
    Fb = sp.Matrix(rng.integers(-2, 3, (N, N)))
    while Fb.det() == 0: Fb = sp.Matrix(rng.integers(-2, 3, (N, N)))
    Esp = [sp.Matrix(E) for E in Es]
    Z = [sp.zeros(N, 0)]; dims = [0]
    for t in range(1, N+1):
        cols = sp.Matrix.hstack(*[E * Fb[:, :t] for E in Esp]); dims.append(cols.rank()); Z.append(cols)
    bound = sum(N - max(t, dims[t]) for t in range(1, N))
    # rows: adapted basis of decreasing chain Z_t^perp; t(b) assignment: bottom rows deepest
    rows = []; tb = []
    for t in range(N-1, -1, -1):
        if t == 0: perp = sp.eye(N)
        else:
            ns = Z[t].T.nullspace(); perp = sp.Matrix.hstack(*ns) if ns else sp.zeros(N, 0)
        for j in range(perp.shape[1]):
            v = perp[:, j]
            cand = sp.Matrix.hstack(*(rows + [v])) if rows else v
            if cand.rank() > len(rows): rows.append(v); tb.append(t)
            if len(rows) == N: break
        if len(rows) == N: break
    rows = rows[::-1]; tb = tb[::-1]            # row 1 = shallowest ... row N = deepest
    X = sp.Matrix.hstack(*rows).T; Yinv = Fb
    Et = [X * E * Yinv for E in Esp]
    Cs = [rng.integers(-2, 3, (N, N)) for _ in range(r)]
    # transformed C: Y C X^{-1}
    Ct = [Yinv.inv() * sp.Matrix(C) * X.inv() for C in Cs]
    # Theta[(b,j),(k,b')] = sum_l Et_l[b,b'] Ct_l[j,k]
    tt = [min([c for c in range(N) if any(E[b, c] != 0 for E in Et)] + [N]) for b in range(N)]
    Om = [(b, j) for b in range(N) for j in range(N) if j < b and j < tt[b]]
    idx = {o: i for i, o in enumerate(Om)}
    T = sp.zeros(len(Om), len(Om))
    for (b, j) in Om:
        for (k, bp) in Om:
            T[idx[(b, j)], idx[(k, bp)]] = sum(Et[l][b, bp] * Ct[l][j, k] for l in range(r))
    nil = (T**max(1, len(Om))).is_zero_matrix if Om else True
    S = Smat(N); K = sum(np.kron(E, C) for E, C in zip(Es, Cs))
    rk = qrank(sp.Matrix(S - K))
    good = nil and len(Om) >= bound and rk >= bound
    tot += 1; ok += good
    print(f"N={N} m={m} r={r} dims={dims} bound={bound} |Om|={len(Om)} nilpotent={nil} rank={rk} {'OK' if good else 'FAIL'}", flush=True)
print("flag lemma:", ok, "/", tot)
