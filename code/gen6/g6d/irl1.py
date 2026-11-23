"""Iteratively reweighted L1 search for sparse per-point errors (commutative class), exact verification.
usage: python3 irl1.py l d [iters] [mode: total|max] [seeds]
Upper bounds only (constructions); verified exactly over Q via milp.exact_feasible."""
import sys, numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix, csr_matrix, hstack, vstack, identity
import milp

def run(l, d, iters=25, mode="total", seed=0):
    S = milp.setup(l, d); taus, lam, V, idx, chi, edges = S
    rng = np.random.default_rng(seed)
    nv = len(V) * d; ne = len(edges)
    lamf = {k: [float(x) for x in lam[k]] for k in lam}
    A = lil_matrix((ne, nv)); c0 = np.zeros(ne)
    for e, (z, w, k) in enumerate(edges):
        for j in range(d):
            A[e, idx[z] * d + j] += lamf[k][j + 1]; A[e, idx[w] * d + j] -= lamf[k][j + 1]
        c0[e] = lamf[k][0] * (chi[z] - chi[w])
    A = csr_matrix(A)
    scale = np.array([np.abs(A[e]).sum() + abs(c0[e]) for e in range(ne)])
    W = rng.uniform(0.5, 1.5, ne)
    best = None
    for it in range(iters):
        # vars: x (nv, free), s (ne >= 0): -s <= (A x + c0)/scale <= s ; min sum W s
        As = csr_matrix(A.multiply(1 / scale[:, None]))
        I = identity(ne, format="csr")
        Aub = vstack([hstack([As, -I]), hstack([-As, -I])])
        bub = np.concatenate([-c0 / scale, c0 / scale])
        cost = np.concatenate([np.zeros(nv), W])
        res = linprog(cost, A_ub=Aub, b_ub=bub, bounds=[(None, None)] * nv + [(0, None)] * ne, method="highs")
        if res.status != 0: break
        r = np.abs(A @ res.x[:nv] + c0) / scale
        bad = {e for e in range(ne) if r[e] > 1e-7}
        per = [0] * l
        for e in bad: per[edges[e][2]] += 1
        val = len(bad) if mode == "total" else (max(per), len(bad))
        if best is None or val < best[0]:
            if milp.exact_feasible(l, d, bad, S):
                best = (val, bad, per)
        eps = 1e-3 if it < iters // 2 else 1e-5
        W = 1.0 / (r + eps)
        if mode == "max":  # penalise heavy directions
            dirw = np.array([per[edges[e][2]] for e in range(ne)], float)
            W *= (1 + dirw / max(1, max(per)))
        W *= rng.uniform(0.9, 1.1, ne)
    return best

if __name__ == "__main__":
    l, d = int(sys.argv[1]), int(sys.argv[2])
    iters = int(sys.argv[3]) if len(sys.argv) > 3 else 25
    mode = sys.argv[4] if len(sys.argv) > 4 else "total"
    seeds = int(sys.argv[5]) if len(sys.argv) > 5 else 3
    bests = [run(l, d, iters, mode, s) for s in range(seeds)]
    bests = [b for b in bests if b]
    b = min(bests, key=lambda x: x[0])
    tot = len(b[1]); conj = (l - d) * 2 ** (l - d - 1)
    print(f"l={l} d={d} mode={mode} best_total={tot} per_dir={b[2]} max_dir={max(b[2])} "
          f"frac_total={tot/(l*2**(l-1)):.4f} max_rho={max(b[2])/2**l:.4f} conj(l-d)2^(l-d-1)={conj} 2^-(d+1)={2**-(d+1):.4f}", flush=True)
