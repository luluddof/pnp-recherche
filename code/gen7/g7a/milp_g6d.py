"""G6d: per-point-error deformation, commutative (Hadamard-diagonal) class.
J(t) = g_t(X_1..X_l) <-> family of polys g_z(t), z in {+-1}^l, deg <= d, g_z(0) = chi(z).
Point tau_k = k is k-nice up to error E_k; rank E_k = #bad edges in direction k,
edge (z, z^k) bad iff g_z(tau_k) != g_{z^k}(tau_k).
MILP: minimise total bad edges (or max per direction). Result verified exactly over Q.
usage: python3 milp.py l d [total|max] [timelimit]
"""
import sys, itertools, numpy as np
from fractions import Fraction as F
from scipy.optimize import milp, LinearConstraint, Bounds
import flint
Q = lambda x: flint.fmpq(x.numerator, x.denominator)

def lagr(nodes, x):
    out = []
    for i, a in enumerate(nodes):
        c = F(1)
        for j, b in enumerate(nodes):
            if j != i: c *= F(x - b, a - b)
        out.append(c)
    return out

def setup(l, d, taus=None):
    taus = taus or list(range(1, l + 1))
    nodes = [0] + taus[:d]
    lam = {k: lagr(nodes, taus[k]) for k in range(l)}   # g(tau_k) = lam0*chi + sum_j lam_j y_j
    V = list(itertools.product([1, -1], repeat=l))
    idx = {z: i for i, z in enumerate(V)}
    chi = {z: int(np.prod(z)) for z in V}
    edges = []
    for z in V:
        for k in range(l):
            if z[k] == 1:
                w = list(z); w[k] = -1; edges.append((z, tuple(w), k))
    return taus, lam, V, idx, chi, edges

def exact_feasible(l, d, bad, S=None):
    """bad: set of edge indices. Exact rational check that good-edge equalities are solvable."""
    taus, lam, V, idx, chi, edges = S or setup(l, d)
    nv = len(V) * d
    rows, rhs = [], []
    for e, (z, w, k) in enumerate(edges):
        if e in bad: continue
        r = [F(0)] * nv
        for j in range(d):
            r[idx[z] * d + j] += lam[k][j + 1]
            r[idx[w] * d + j] -= lam[k][j + 1]
        rows.append(r); rhs.append(-lam[k][0] * (chi[z] - chi[w]))
    if not rows: return True
    A = flint.fmpq_mat(len(rows), nv, [Q(x) for r in rows for x in r])
    Ab = flint.fmpq_mat(len(rows), nv + 1, [Q(x) for r, b in zip(rows, rhs) for x in r + [b]])
    return A.rank() == Ab.rank()

def solve(l, d, mode="total", tl=300, Y=30.0, perfect=0):
    S = setup(l, d); taus, lam, V, idx, chi, edges = S
    nv = len(V) * d; ne = len(edges); nl = l
    lamf = {k: [float(x) for x in lam[k]] for k in lam}
    M = {k: 2 * abs(lamf[k][0]) + 2 * Y * sum(abs(x) for x in lamf[k][1:]) + 1 for k in range(l)}
    nvar = nv + ne + (1 if mode == "max" else 0)
    A = []; lo = []; hi = []
    for e, (z, w, k) in enumerate(edges):
        r = np.zeros(nvar)
        for j in range(d):
            r[idx[z] * d + j] += lamf[k][j + 1]; r[idx[w] * d + j] -= lamf[k][j + 1]
        c0 = lamf[k][0] * (chi[z] - chi[w])
        r1 = r.copy(); r1[nv + e] = -M[k]; A.append(r1); lo.append(-np.inf); hi.append(-c0)
        r2 = r.copy(); r2[nv + e] = M[k]; A.append(r2); lo.append(-c0); hi.append(np.inf)
    c = np.zeros(nvar)
    if mode == "max":
        for k in range(l):
            r = np.zeros(nvar); r[-1] = -1
            for e, (z, w, kk) in enumerate(edges):
                if kk == k: r[nv + e] = 1
            A.append(r); lo.append(-np.inf); hi.append(0)
        c[-1] = 1
    else:
        c[nv:nv + ne] = 1
    integ = np.zeros(nvar); integ[nv:nv + ne] = 1
    lb = np.full(nvar, -Y); ub = np.full(nvar, Y); lb[nv:] = 0; ub[nv:nv + ne] = 1
    for e, (z, w, k) in enumerate(edges):
        if k < perfect: ub[nv + e] = 0
    if mode == "max": ub[-1] = np.inf
    res = milp(c, constraints=LinearConstraint(np.array(A), lo, hi), integrality=integ,
               bounds=Bounds(lb, ub), options={"time_limit": tl, "disp": False})
    b = res.x[nv:nv + ne]
    bad = {e for e in range(ne) if b[e] > 0.5}
    per = [sum(1 for e in bad if edges[e][2] == k) for k in range(l)]
    ok = exact_feasible(l, d, bad, S)
    return res, bad, per, ok

def subcube_lb_total(l, d):
    from math import comb
    best = 0
    for s in range(d + 1, l + 1):
        best = max(best, -(-l * (s - d) * 2 ** (l - s) // s))
    return best

if __name__ == "__main__":
    l, d = int(sys.argv[1]), int(sys.argv[2])
    mode = sys.argv[3] if len(sys.argv) > 3 else "total"
    tl = float(sys.argv[4]) if len(sys.argv) > 4 else 300
    pf = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    res, bad, per, ok = solve(l, d, mode, tl, perfect=pf)
    tot = len(bad)
    print(f"l={l} d={d} mode={mode} status={res.status} obj={res.fun:.3f} dual_bound={getattr(res,'mip_dual_bound',None)} "
          f"total_bad={tot} per_dir={per} exact_ok={ok} frac={tot/(l*2**(l-1)):.4f} subcubeLB_total={subcube_lb_total(l,d)}")

def show(l, d, mode="total", tl=60):
    S = setup(l, d); edges = S[5]
    res, bad, per, ok = solve(l, d, mode, tl)
    print(f"l={l} d={d} {mode} total={len(bad)} per={per} ok={ok} status={res.status}")
    for e in sorted(bad, key=lambda e: edges[e][2]):
        z, w, k = edges[e]
        print("  dir", k + 1, "z(other coords)=", "".join('+' if z[j] == 1 else '-' for j in range(l) if j != k))
