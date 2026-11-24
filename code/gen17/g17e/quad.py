"""G18d: dual (quadrature) form. Nodes = edges (z,k) (z_k=0). W_{D'} basis chi_S(z) pi_S(t) t^i, |S|+i<=D'=l-d-1.
B feasible <=> lambda (coefficient of (S=empty,i=0)) in span{ev_e: e in B}.
Heuristic: iteratively reweighted L1 on the weights w (float LP), then EXACT verification mod p:
 (a) lambda in rowspan(M_B) mod p, (b) primal system on good edges consistent mod p (independent check, small l)."""
import sys, itertools, random, time
import numpy as np, flint
from math import comb
from scipy.optimize import linprog
import scipy.sparse as sp
P = 2**31 - 1
def taus_for(l, kind="int"):
    if kind == "int": return [k + 1 for k in range(l)]
    if kind == "sym": return [(k // 2 + 1) * (1 if k % 2 == 0 else -1) for k in range(l)]
    rng = random.Random(5); return [rng.randrange(2, 50) * rng.choice([1, -1]) for _ in range(l)]
def basis(l, Dp): return [(S, i) for s in range(Dp + 1) for S in itertools.combinations(range(l), s) for i in range(Dp - s + 1)]
def nodes(l): return [(z, k) for z in range(1 << l) for k in range(l) if not z >> k & 1]
def bound(l, d): return sum(comb(l, j) * max(0, l - d - 2 * j) for j in range(l + 1))
def layered(l, d):
    best = None
    for shift in range(d + 1):
        cuts = sorted(set([0] + list(range(shift if shift > 0 else d + 1, l + 1, d + 1)) + [l + 1])); bnd = cuts[1:-1]
        B = [(z, k) for z, k in nodes(l) if (bin(z).count('1') + 1) in bnd]
        if best is None or len(B) < len(best): best = B
    return best
def mat_float(l, d, tau):
    Dp = l - d - 1; bs = basis(l, Dp); nd = nodes(l)
    Z = np.array([[1 - 2 * (z >> j & 1) for j in range(l)] for z, k in nd], dtype=float)
    K = np.array([k for z, k in nd]); T = np.array([tau[k] for k in K], dtype=float)
    U = Z * (1 - T[:, None] / np.array(tau, dtype=float)[None, :])   # u_j at node (u_k = 0 automatically)
    cols = []
    for S, i in bs:
        v = T ** i
        for j in S: v = v * U[:, j]
        cols.append(v)
    return np.array(cols).T, bs, nd
def row_mod(l, tau, inv, bs, z, k):
    t = tau[k] % P; r = []
    for S, i in bs:
        if k in S: r.append(0); continue
        v = pow(t, i, P)
        for j in S:
            v = v * ((1 - t * inv[j]) % P) % P
            if z >> j & 1: v = (-v) % P
        r.append(v)
    return r
def exact_dual(l, d, tau, B):
    Dp = l - d - 1; bs = basis(l, Dp); inv = [pow(t % P, P - 2, P) for t in tau]
    R = [row_mod(l, tau, inv, bs, z, k) for z, k in B]
    lam = [1 if (len(S) == 0 and i == 0) else 0 for S, i in bs]
    r1 = flint.nmod_mat(R, P).rank() if R else 0
    r2 = flint.nmod_mat(R + [lam], P).rank()
    return r1, r1 == r2
def exact_primal(l, d, tau, B):
    Bs = set(B); n = (1 << l) * d; R = []
    for z, k in nodes(l):
        if (z, k) in Bs: continue
        w = z | 1 << k; r = [0] * (n + 1); tk = tau[k] % P; pw = tk
        for j in range(d):
            r[z * d + j] = pw; r[w * d + j] = pw; pw = pw * tk % P
        r[n] = P - 2; R.append(r)
    if not R: return True
    return flint.nmod_mat([r[:n] for r in R], P).rank() == flint.nmod_mat(R, P).rank()
def per_dir(l, B): return [sum(1 for z, k in B if k == kk) for kk in range(l)]
def irl1(l, d, tau, iters=12, seed=0, mode="total", verbose=False):
    M, bs, nd = mat_float(l, d, tau); ne, nb = M.shape
    sc = np.abs(M).max(axis=0); sc[sc == 0] = 1; Ms = M / sc
    lam = np.zeros(nb); lam[0] = 1.0
    assert bs[0] == ((), 0)
    rng = np.random.default_rng(seed); c = np.ones(ne) + (0.3 * rng.random(ne) if seed else 0)
    best = None
    A = sp.csr_matrix(np.hstack([Ms.T, -Ms.T]))
    for it in range(iters):
        res = linprog(np.concatenate([c, c]), A_eq=A, b_eq=lam, bounds=(0, None), method="highs")
        if res.status != 0: break
        w = res.x[:ne] - res.x[ne:]; mx = np.abs(w).max()
        supp = [nd[e] for e in range(ne) if abs(w[e]) > 1e-7 * mx]
        if best is None or len(supp) < len(best): best = supp
        if verbose: print("   it", it, len(supp), flush=True)
        c = 1.0 / (np.abs(w) / mx + 1e-3)
    return best
if __name__ == "__main__":
    cases = [tuple(map(int, a.split(","))) for a in sys.argv[2:]]
    kind = sys.argv[1]
    for l, d in cases:
        tau = taus_for(l, kind); t0 = time.time()
        L = layered(l, d); rkL, okL = exact_dual(l, d, tau, L)
        best = None
        for seed in range(4):
            B = irl1(l, d, tau, seed=seed)
            if B is None: continue
            rk, ok = exact_dual(l, d, tau, B)
            pr = exact_primal(l, d, tau, B) if l <= 8 else None
            if ok and (best is None or len(B) < len(best[0])): best = (B, rk, pr)
            if not ok: print("   seed", seed, "support", len(B), "NOT exact-feasible")
        print(f"(l,d)=({l},{d}) tau={kind}: bound={bound(l,d)} layered={len(L)} rank_layered={rkL} ({okL}) "
              f"IRL1 best={len(best[0]) if best else None} rank={best[1] if best else None} primal_ok={best[2] if best else None} "
              f"per_dir={per_dir(l,best[0]) if best else None} {time.time()-t0:.0f}s", flush=True)
