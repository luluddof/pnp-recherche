# Large-h version: F^mu_e with c = n (no sketch of sigma'), streamed per monomial and projected by a DENSE random
# matrix S (T x n^3) mod P (faithful w.h.p. when the projected rank of [all cols | rhs] is < T).
import sys, time, pickle
sys.dont_write_bytecode = True
import numpy as np
from math import comb
from appb_exact import partitions
from fcomp import Ops, mm

def proj(S, V, P):
    """S (T x m float64, entries < P < 2^24), V (m x k int64) -> S @ V mod P exactly (m < 2^16.9)."""
    V1 = (V >> 12).astype(np.float64); V0 = (V & 4095).astype(np.float64)
    r1 = np.fmod(S @ V1, P).astype(np.int64); r0 = np.fmod(S @ V0, P).astype(np.int64)
    return ((r1 << 12) + r0) % P

def run(h, D, T, P=16777213, seed=7, wset=None):
    t0 = time.time(); O = Ops(P); Vh = O.basis(h); n = len(Vh); c = n
    mus = [m for w in range(D + 1) for m in sorted(partitions(w)) if wset is None or w in wset]
    G = np.eye(n, dtype=np.int64)            # c = n: identity = no sketch
    m_rows = n * n * c
    rng = np.random.default_rng(seed)
    S = rng.integers(0, P, size=(T, m_rows)).astype(np.float64)
    wm = max(2 * h, h + D)
    BpG = {}; BmS = {}
    for w in range(0, wm + 1):
        if w <= 2 * h: BpG[w] = np.concatenate([O.mat(nu, w - h, h) for nu in Vh], axis=1)   # (pw, n*n) [nu, sigma']
        BmS[w] = np.concatenate([O.mat(nu, h - w, w) for nu in Vh], axis=0)                   # (n*n, pw) [nu, sigma]
    print(f"h={h} D={D} n={n} rows={m_rows} T={T}: b-operators {time.time()-t0:.0f}s", flush=True)
    cols = {}; bad = 0
    for im, mu in enumerate(mus):
        hp = sum(mu); H = h + hp; fp = {}; fm = {}
        for w in range(0, wm + 1):
            if w <= 2 * h:
                X = mm(O.mat(mu, h - w, w), BpG[w], P)                     # (sigma, [nu, sigma'])
                fp[w - 2 * h] = X.reshape(n, n, n).transpose(1, 0, 2).reshape(-1)
            Y = mm(BmS[w], O.mat(mu, w - h, h), P)                          # ([nu, sigma], sigma')
            fm[-w] = Y.reshape(-1)
        def Fe(e, f, lo):
            acc = np.zeros(m_rows, dtype=np.int64)
            for j in range(0, H + 1):
                mm_ = e - j
                if lo <= mm_ <= 0 and mm_ in f:
                    acc = (acc + ((-1) ** j * comb(H, j)) % P * f[mm_]) % P
            return acc
        es = list(range(-2 * h, h + hp + 1))
        Fs = [Fe(e, fp, -2 * h) if e <= 0 else Fe(e, fm, -wm) for e in es]
        if not np.array_equal(Fe(0, fp, -2 * h), Fe(0, fm, -wm)): bad += 1
        Pm = proj(S, np.stack(Fs, axis=1), P)
        for i, e in enumerate(es): cols[(mu, e)] = Pm[:, i].copy()
        if im % 20 == 0: print(f"  mu {im}/{len(mus)} {time.time()-t0:.0f}s", flush=True)
    rhs = []
    for k0 in range(n):
        v = np.zeros((n, n, n), dtype=np.int64); v[k0] = G
        rhs.append(v.reshape(-1))
    rhsP = proj(S, np.stack(rhs, axis=1), P)
    rhs = [rhsP[:, i].copy() for i in range(n)]
    print(f"h={h}: done {len(cols)} cols, e=0 two-route failures {bad}, {time.time()-t0:.0f}s", flush=True)
    return dict(h=h, D=D, c=c, P=P, mus=mus, cols=cols, rhs=rhs, n=n, T=T)

if __name__ == "__main__":
    h, D, T = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    d = run(h, D, T)
    from compress import compress
    out = compress(d)
    print(f"r_all = {out['r_all']}  (T = {T}; faithful w.h.p. iff r_all < T)", flush=True)
    pickle.dump(out, open(f"Fbig_h{h}_D{D}_T{T}_cmp.pkl", "wb"))
