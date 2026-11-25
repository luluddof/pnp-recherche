# G27-T1: k_max(h, D) for the degree-1 ansatz (P1's F-route), exact rows (c = n, no sketch, no projection), mod P.
# usage: python3 dscan.py h D [P] [Wlist]   -> prints, per W, k_max(W) = 4h-1-W-minW' and ranks.
import sys, time, pickle
sys.dont_write_bytecode = True
import numpy as np, flint
from math import comb
from appb_exact import partitions, apply_mono
from fcomp import Ops, mm

def build(h, D, P):
    t0 = time.time(); O = Ops(P); Vh = O.basis(h); n = len(Vh)
    mus = [m for w in range(D + 1) for m in sorted(partitions(w))]
    m_rows = n ** 3; wm = max(2 * h, h + D)
    BpG = {}; BmS = {}
    for w in range(0, wm + 1):
        if w <= 2 * h: BpG[w] = np.concatenate([O.mat(nu, w - h, h) for nu in Vh], axis=1)
        BmS[w] = np.concatenate([O.mat(nu, h - w, w) for nu in Vh], axis=0)
    keys = []; blocks = []; bad = 0
    for im, mu in enumerate(mus):
        hp = sum(mu); H = h + hp; fp = {}; fm = {}
        for w in range(0, wm + 1):
            if w <= 2 * h:
                X = mm(O.mat(mu, h - w, w), BpG[w], P)
                fp[w - 2 * h] = X.reshape(n, n, n).transpose(1, 0, 2).reshape(-1)
            fm[-w] = mm(BmS[w], O.mat(mu, w - h, h), P).reshape(-1)
        def Fe(e, f, lo):
            acc = np.zeros(m_rows, dtype=np.int64)
            for j in range(0, H + 1):
                q = e - j
                if lo <= q <= 0 and q in f: acc = (acc + ((-1) ** j * comb(H, j)) % P * f[q]) % P
            return acc
        es = list(range(-2 * h, h + hp + 1))
        if not np.array_equal(Fe(0, fp, -2 * h), Fe(0, fm, -wm)): bad += 1
        for e in es:
            keys.append((hp, e)); blocks.append(Fe(e, fp, -2 * h) if e <= 0 else Fe(e, fm, -wm))
        if hp >= 10: apply_mono.cache_clear()
    rhs = []
    for k0 in range(n):
        v = np.zeros((n, n, n), dtype=np.int64); v[k0] = np.eye(n, dtype=np.int64); rhs.append(v.reshape(-1))
    A = np.stack(blocks, axis=1) % P
    print(f"h={h} D={D} P={P}: n={n} rows={m_rows} cols={A.shape[1]} two-route failures {bad}  {time.time()-t0:.0f}s", flush=True)
    return n, keys, A, np.stack(rhs, axis=1)

def rank(A, P):
    if A.shape[1] == 0: return 0
    return flint.nmod_mat(A.shape[0], A.shape[1], A.ravel().tolist(), P).rank()

def compress_rows(A, R, P):
    """keep a row subset on which [A|R] keeps its rank (exact: row-basis of [A|R])."""
    M = np.concatenate([A, R], axis=1)
    MT = flint.nmod_mat(M.shape[1], M.shape[0], M.T.ravel().tolist(), P)
    Rr, r = MT.rref()
    piv = []
    for i in range(r):
        for j in range(M.shape[0]):
            if int(Rr[i, j]): piv.append(j); break
    piv = np.array(piv)
    return A[piv], R[piv], r

def scan(h, D, P=16777213, Ws=None):
    n, keys, A, R = build(h, D, P)
    t = time.time()
    # reduce columns first (many columns): row compression of [A|R]
    A, R, r_all = compress_rows(A, R, P)
    print(f"  r_all = {r_all} (of n^3={n**3}); compressed {time.time()-t:.0f}s", flush=True)
    hp = np.array([k[0] for k in keys]); ee = np.array([k[1] for k in keys])
    def feas(W, Wp):
        sel = (ee <= Wp - 2 * h) | (ee >= hp + h - W)
        As = A[:, sel]; ra = rank(As, P); rb = rank(np.concatenate([As, R], axis=1), P)
        return ra == rb, ra, rb
    out = {}
    Ws = Ws if Ws is not None else list(range(0, 4 * h))
    for W in Ws:
        lo, hi = -1, 4 * h - 1 - W      # W' = -1 : kappa = 0
        if not feas(W, hi)[0]: print(f"  W={W}: threshold infeasible?!"); continue
        while lo < hi:
            mid = (lo + hi) // 2
            if feas(W, mid)[0]: hi = mid
            else: lo = mid + 1
        f = feas(W, lo); g = feas(W, lo - 1) if lo > -1 else (False, 0, 0)
        k = 4 * h - 1 - W - lo; out[W] = (lo, k, f[1], g[2] - g[1])
        print(f"  W={W:>3}: min W'={lo:>3} k={k:>3}  rank(allowed)={f[1]} (r_all {r_all}); at W'-1 deficiency {g[2]-g[1]}  {time.time()-t:.0f}s", flush=True)
    return dict(h=h, D=D, P=P, r_all=r_all, n=n, res=out)

if __name__ == "__main__":
    h, D = int(sys.argv[1]), int(sys.argv[2]); P = int(sys.argv[3]) if len(sys.argv) > 3 else 16777213
    Ws = [int(a) for a in sys.argv[4].split(',')] if len(sys.argv) > 4 else None
    res = scan(h, D, P, Ws)
    pickle.dump(res, open(f"dscan_h{h}_D{D}_P{P}.pkl", "wb"))
