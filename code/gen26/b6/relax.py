# Test A/B: "best perturbation inside App B's factor spaces".
# U_S  = span{ g_sigma : sigma in S },  g_sigma(u,v)  = [t^mu(u)] b_v[h-|sigma|] t^sigma      (row side of R)
# U'_S' = span{ g'_sigma: sigma in S'}, g'_sigma(y,v) = [t^sigma] b_v[|sigma|-h] t^mu(y)      (col side of C)
# Question: is X = R + C with R in Q^{N^2}_(xy) (x) U_S, C in Q^{N^2}_(xu) (x) U'_S' ?   (x-side completely free)
# Reduction: K = annihilator of U_S (in Q^{N^2}_(uv)); need C with (I (x) K) C = (I (x) K) X; decouples over x:
#   M[(j,y),(m,u)] = (G'_m K_j^T)[y,u],  RHS block j = K_j (columns x).  Feasible iff rank M = rank [M | RHS].
# All ranks mod a large prime (two primes for safety).
import sys, math, time, json, itertools
sys.dont_write_bytecode = True
import flint
from appb_exact import partitions, apply_mono, zmu
PR = [33554393, 16777213]
import numpy as np

def setup(h, N=None):
    parts = sorted(partitions(h)); N = N or len(parts); mus = parts[:N]
    bs = [math.factorial(h) // zmu(nu) for nu in mus]
    return N, mus, bs

def g_rows(h, N, mus, bs, wmax):
    """dict sigma -> vector over (u,v) (row side, U) for |sigma| <= wmax"""
    idx = {m: i for i, m in enumerate(mus)}
    out = {}
    for w in range(0, wmax + 1):
        for sig in partitions(w):
            vec = [0] * (N * N)
            for v, nu in enumerate(mus):
                for s2, c in apply_mono(nu, h - w, sig):
                    u = idx.get(s2)
                    if u is not None: vec[u * N + v] += bs[v] * c
            out[sig] = vec
    return out

def gp_rows(h, N, mus, bs, wmax):
    """dict sigma -> vector over (y,v) (col side, U')"""
    out = {}
    for w in range(0, wmax + 1):
        for sig in partitions(w): out[sig] = [0] * (N * N)
        for v, nu in enumerate(mus):
            for y, muy in enumerate(mus):
                for s2, c in apply_mono(nu, w - h, muy):
                    out[s2][y * N + v] += bs[v] * c
    return out

def basis_mod(vecs, n, p):
    if not vecs: return []
    M = flint.nmod_mat(len(vecs), n, [x % p for v in vecs for x in v], p)
    R = M.rref()[0]; r = M.rank()
    return [[int(R[i, j]) for j in range(n)] for i in range(r)]

def annihilator(vecs, n, p):
    if not vecs: return [[int(i == j) for j in range(n)] for i in range(n)]
    M = flint.nmod_mat(len(vecs), n, [x % p for v in vecs for x in v], p)
    X, nul = M.nullspace()
    return [[int(X[i, j]) for i in range(n)] for j in range(nul)]

def feasible(N, Uvecs, Upvecs, p):
    n = N * N
    K = annihilator(Uvecs, n, p)
    if not K: return True, n, len(basis_mod(Upvecs, n, p))
    Bp = basis_mod(Upvecs, n, p)
    dU = n - len(K); dUp = len(Bp)
    rows = len(K) * N; cols = dUp * N
    Ka = np.array(K, dtype=np.int64).reshape(len(K), N, N) % p      # [j,u,v]
    Ga = np.array(Bp, dtype=np.int64).reshape(dUp, N, N) % p if dUp else np.zeros((0, N, N), dtype=np.int64)
    Mt = np.einsum('myv,juv->jymu', Ga, Ka) % p                      # [j,y,m,u]
    Mm = Mt.reshape(rows, cols)
    Rhs = Ka.reshape(rows, N)                                        # block j = K_j[y][x]  (u->y, v->x)
    ent = np.concatenate([Mm, Rhs], axis=1).tolist()
    A = flint.nmod_mat(rows, cols + N, [e for row in ent for e in row], p)
    if cols == 0: return A.rank() == 0, dU, dUp
    M = flint.nmod_mat(rows, cols, [ent[r][c] for r in range(rows) for c in range(cols)], p)
    return M.rank() == A.rank(), dU, dUp

if __name__ == "__main__":
    h = int(sys.argv[1]); N = int(sys.argv[2]) if len(sys.argv) > 2 else None
    t0 = time.time()
    N, mus, bs = setup(h, N)
    wmax = 4 * h
    G = g_rows(h, N, mus, bs, wmax); Gp = gp_rows(h, N, mus, bs, wmax)
    print(f"h={h} N={N}: vectors built ({time.time()-t0:.1f}s)", flush=True)
    p = PR[0]; n = N * N
    # dims by weight
    dimU = [len(basis_mod([G[s] for s in G if sum(s) <= W], n, p)) for W in range(wmax + 1)]
    dimUp = [len(basis_mod([Gp[s] for s in Gp if sum(s) <= W], n, p)) for W in range(wmax + 1)]
    print("dim U_W :", dimU); print("dim U'_W:", dimUp, flush=True)
    res = []
    for W in range(-1, wmax + 1):
        Uv = [G[s] for s in G if sum(s) <= W]
        if W >= 0 and dimU[W] == n: res.append((W, -1, n, 0)); print(W, "U full"); break
        best = None
        lo, hi = -1, prevWp if 'prevWp' in dir() else wmax
        def test(Wp):
            Upv = [Gp[s] for s in Gp if sum(s) <= Wp]
            return feasible(N, Uv, Upv, p)
        while lo < hi:
            mid = (lo + hi) // 2
            if test(mid)[0]: hi = mid
            else: lo = mid + 1
        ok, dU, dUp = test(lo); ok2 = feasible(N, Uv, [Gp[s] for s in Gp if sum(s) <= lo], PR[1])[0]
        okm = test(lo - 1)[0] if lo > -1 else False
        best = (W, lo, dU, dUp, ok2 and ok and not okm); prevWp = lo
        thr = 4 * h - 1 - W
        print(f"W={W:>3}: min W'={best[1]:>3} (threshold W'={thr}), dims (rr,rc) <= ({best[2]},{best[3]})  checks ok={best[4]}", flush=True)
        res.append(best[:4])
    json.dump({"h": h, "N": N, "dimU": dimU, "dimUp": dimUp, "res": res}, open(f"relax_h{h}_N{N}.json", "w"))
    print(f"done {time.time()-t0:.1f}s")
