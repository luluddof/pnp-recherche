"""G5a experiment: crossing expansions over product vs non-product partition systems.

Object (generalises OpenAI's J(t) / G4c lemma):
  V finite index set, partitions pi_1..pi_D of V  (phi: V -> prod [q_i]),
  s(x,y) = #{i : pi_i(x) != pi_i(y)}  (Hamming metric of the embedding),
  crossing C = permutation of V with s(x, Cx) = D (moves every class).
  Want J(t) in Q[t]^{VxV}, deg_t J <= d, J(0) = C, J(t_i) = K_i + E for D points t_i,
  K_i pi_i-block-diagonal, E common error.  Gain = d/D (Taylor order h = n d / D).

Exact characterisation: such J exists with error E  <=>  E[x,y] = C[x,y] for every pair with
s(x,y) > d (near pairs free: interpolate s(x,y)+1 <= d+1 values).  So
  r*(d) = min rank of a completion of C on the far pattern  {(x,y): s(x,y) > d}.
Lower bound: alpha(H_d), H_d: a~a' iff s(a,Ca')<=d or s(a',Ca)<=d   (then E[A,CA] = I).
Upper bound for C = translation on F_2^k (linear code G): best translation-invariant E,
  = min Fourier support of f with f(w)=1, f(z)=0 for z != w with wt(zG) > d  (MILP).
"""
import itertools, sys, json
import numpy as np
import networkx as nx
from scipy.optimize import milp, LinearConstraint, Bounds

rng = np.random.default_rng(1)


def words(k):
    return np.array(list(itertools.product([0, 1], repeat=k)), dtype=np.int64)


def alpha_lower(phi, Cperm, d, tl=40):
    """independence number of H_d via MILP (returns (value, proven_optimal))."""
    K = len(phi)
    s = (phi[:, None, :] != phi[None, :, :]).sum(-1)
    sC = s[:, Cperm]
    adj = (sC <= d) | (sC.T <= d)
    np.fill_diagonal(adj, False)
    E = np.argwhere(np.triu(adj))
    if len(E) == 0:
        return K, True
    A = np.zeros((len(E), K)); A[np.arange(len(E)), E[:, 0]] = 1; A[np.arange(len(E)), E[:, 1]] = 1
    res = milp(-np.ones(K), constraints=[LinearConstraint(A, -np.inf, 1)], integrality=np.ones(K),
               bounds=Bounds(0, 1), options={"time_limit": tl})
    # report the dual bound if not optimal is not needed: a feasible independent set is a valid LB
    return int(round(-res.fun)), res.status == 0


def ti_rank_binary(G, w, d, M=64.0):
    """min #nonzero Fourier coeffs of f: F_2^k -> R with f(w)=1, f(z)=0 for far z != w."""
    k, D = G.shape
    Z = words(k)
    wts = (Z @ G % 2).sum(1)
    widx = int(np.where((Z == w).all(1))[0][0])
    far = [i for i in range(len(Z)) if i != widx and wts[i] > d]
    rows = [widx] + far
    H = (-1.0) ** ((Z @ Z.T) % 2)  # H[z,s] = chi_s(z)
    N = len(Z)
    # variables: fhat (N, free in [-M,M]), y (N, binary)
    c = np.concatenate([np.zeros(N), np.ones(N)])
    A_eq = np.hstack([H[rows, :], np.zeros((len(rows), N))])
    b_eq = np.zeros(len(rows)); b_eq[0] = 1.0
    # |fhat_s| <= M y_s
    I = np.eye(N)
    A1 = np.hstack([I, -M * I]); A2 = np.hstack([-I, -M * I])
    cons = [LinearConstraint(A_eq, b_eq, b_eq),
            LinearConstraint(np.vstack([A1, A2]), -np.inf, 0.0)]
    integ = np.concatenate([np.zeros(N), np.ones(N)])
    bnds = Bounds(np.concatenate([-M * np.ones(N), np.zeros(N)]),
                  np.concatenate([M * np.ones(N), np.ones(N)]))
    res = milp(c, constraints=cons, integrality=integ, bounds=bnds,
               options={"time_limit": 60})
    if res.x is None:
        return None, False
    return int(round(res.x[N:].sum())), res.status == 0


def code_structure(G):
    """V = F_2^k, phi(x) = xG, C = translation by w with wG = all ones."""
    k, D = G.shape
    Z = words(k)
    phi = Z @ G % 2
    w = None
    for z in Z:
        if ((z @ G) % 2).all():
            w = z; break
    assert w is not None, "code must contain the all-ones word"
    idx = {tuple(z): i for i, z in enumerate(Z)}
    Cperm = np.array([idx[tuple((z + w) % 2)] for z in Z])
    return phi, Cperm, w


def affine_plane(p):
    """V = F_p^2, partitions = parallel classes except the one of w (C = x + w)."""
    pts = [(a, b) for a in range(p) for b in range(p)]
    dirs = [(1, m) for m in range(p)] + [(0, 1)]
    w = (0, 1)  # translation along direction (0,1); drop that class
    dirs = [dd for dd in dirs if dd != (0, 1)]
    # class of point x for direction (u,v): the line value v*a - u*b mod p
    phi = np.array([[(dv * a - du * b) % p for (du, dv) in dirs] for (a, b) in pts])
    idx = {pt: i for i, pt in enumerate(pts)}
    Cperm = np.array([idx[((a + w[0]) % p, (b + w[1]) % p)] for (a, b) in pts])
    return phi, Cperm


def rm1(m):
    """Reed-Muller RM(1,m): [2^m, m+1], weights 0, 2^{m-1}, 2^m."""
    pts = words(m).T  # m x 2^m
    return np.vstack([np.ones((1, 2 ** m), dtype=np.int64), pts])


def random_code(k, D):
    while True:
        G = rng.integers(0, 2, size=(k, D))
        G[0, :] = 1  # contains all-ones
        if np.linalg.matrix_rank(G) == k and len({tuple(c) for c in G.T}) == D:
            # check rank over F_2
            Z = words(k)
            if len({tuple(r) for r in (Z @ G % 2)}) == 2 ** k:
                return G


def run(name, G=None, phi=None, Cperm=None):
    w = None
    if G is not None:
        phi, Cperm, w = code_structure(G)
    K, D = phi.shape
    out = []
    for d in range(0, D):
        a, aopt = alpha_lower(phi, Cperm, d)
        t, topt = ti_rank_binary(G, w, d) if G is not None else (None, None)
        out.append((d, a, aopt, t, topt))
        if t == 1 or (G is None and a == 1):
            break
    lines = [f"\n== {name}: K={K}, D={D}", "  d  gain=d/D  alpha_LB(*=proven max)  TI_rank_UB(*=proven TI-opt)  [fraction of K]"]
    for d, a, aopt, t, topt in out:
        ts = f"{t:4d}{'*' if topt else ' '} ({t/K:.3f})" if t is not None else "   -"
        lines.append(f"  {d:2d}  {d/D:.3f}    {a:4d}{'*' if aopt else ' '} ({a/K:.3f})   {ts}")
    print("\n".join(lines)); sys.stdout.flush()
    return dict(name=name, K=K, D=D, rows=out)


CASES = {
  "cube4": lambda: run("cube Q_4 (product)", G=np.eye(4, dtype=np.int64)),
  "cube5": lambda: run("cube Q_5 (product)", G=np.eye(5, dtype=np.int64)),
  "cube6": lambda: run("cube Q_6 (product)", G=np.eye(6, dtype=np.int64)),
  "rm13": lambda: run("RM(1,3) [8,4]", G=rm1(3)),
  "rm14": lambda: run("RM(1,4) [16,5]", G=rm1(4)),
  "rm15": lambda: run("RM(1,5) [32,6]", G=rm1(5)),
  "rnd8_4": lambda: run("random [8,4] code w/ all-ones", G=random_code(4, 8)),
  "rnd10_5": lambda: run("random [10,5] code w/ all-ones", G=random_code(5, 10)),
  "rnd15_5": lambda: run("random [15,5] code w/ all-ones", G=random_code(5, 15)),
  "rnd12_6": lambda: run("random [12,6] code w/ all-ones", G=random_code(6, 12)),
  "ag3": lambda: run("affine plane AG(2,3) (2-design)", *([None] + list(affine_plane(3)))),
  "ag5": lambda: run("affine plane AG(2,5) (2-design)", *([None] + list(affine_plane(5)))),
  "ag7": lambda: run("affine plane AG(2,7) (2-design)", *([None] + list(affine_plane(7)))),
}

if __name__ == "__main__":
    r = CASES[sys.argv[1]]()
    json.dump(r, open(f"res_{sys.argv[1]}.json", "w"))
