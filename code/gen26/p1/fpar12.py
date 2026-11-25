# Parity-split + chunked version of fbig.py.  The automorphism t_k -> -t_k (alpha -> -alpha of the Heisenberg VOA) gives
# <sigma| t^mu[r] t^nu[-r] |sigma'> = 0 unless l(mu)+l(nu)+l(sigma)+l(sigma') is even (l = number of parts).
# Sector pi: monomials mu with l(mu) = pi mod 2, rows (nu,sigma,sigma') with l(nu)+l(sigma)+l(sigma') = pi, rhs l = e_nu0^*
# with l(nu0) = pi.  The full problem is feasible iff both sectors are.
# usage: python3 fpar.py h D pi T chunk nchunks   -> part file;   python3 fpar.py merge h D pi T nchunks
import sys, time, pickle
sys.dont_write_bytecode = True
import numpy as np
from math import comb
from appb_exact import partitions
from fcomp import Ops, mm
def proj(S, V, P):
    V1 = (V >> 12).astype(np.float64); V0 = (V & 4095).astype(np.float64); out = []
    for i in range(0, S.shape[0], 256):
        Sb = S[i:i + 256].astype(np.float64)
        r1 = np.fmod(Sb @ V1, P).astype(np.int64); r0 = np.fmod(Sb @ V0, P).astype(np.int64)
        out.append(((r1 << 12) + r0) % P)
    return np.concatenate(out, axis=0)
P = 16777213

def setup(h, D, pi):
    O = Ops(P); Vh = O.basis(h); n = len(Vh)
    lp = np.array([len(s) % 2 for s in Vh])
    par = (lp[:, None, None] + lp[None, :, None] + lp[None, None, :]) % 2      # [nu, sigma, sigma']
    rows = np.nonzero(par.reshape(-1) == pi)[0]
    mus = [m for w in range(D + 1) for m in sorted(partitions(w)) if len(m) % 2 == pi]
    return O, Vh, n, rows, mus, lp

def part(h, D, pi, T, chunk, nchunks, seed=11):
    t0 = time.time(); O, Vh, n, rows, mus, lp = setup(h, D, pi)
    mine = mus[chunk::nchunks] if chunk >= 0 else []
    wm = max(2 * h, h + D); BpG = {}; BmS = {}
    import os
    cache = f"bops_h{h}_D{D}.pkl"
    if chunk == -2:      # partial b-ops for w in [nchunks // 100, nchunks % 100]
        lo, hi = nchunks // 100, nchunks % 100
        for w in range(lo, hi + 1):
            if w <= 2 * h: BpG[w] = np.concatenate([O.mat(nu, w - h, h) for nu in Vh], axis=1)
            BmS[w] = np.concatenate([O.mat(nu, h - w, w) for nu in Vh], axis=0)
        pickle.dump((BpG, BmS), open(f"bops_h{h}_D{D}_w{lo}-{hi}.pkl", "wb")); return
    if chunk == -3:      # combine partial files
        import glob
        for f in glob.glob(f"bops_h{h}_D{D}_w*.pkl"):
            a, b = pickle.load(open(f, "rb")); BpG.update(a); BmS.update(b)
        assert sorted(BmS) == list(range(wm + 1)), sorted(BmS)
        pickle.dump((BpG, BmS), open(cache, "wb")); return
    if os.path.exists(cache):
        BpG, BmS = pickle.load(open(cache, "rb"))
    else:
        for w in range(0, wm + 1):
            if w <= 2 * h: BpG[w] = np.concatenate([O.mat(nu, w - h, h) for nu in Vh], axis=1)
            BmS[w] = np.concatenate([O.mat(nu, h - w, w) for nu in Vh], axis=0)
        pickle.dump((BpG, BmS), open(cache + ".tmp", "wb")); os.replace(cache + ".tmp", cache)
    if chunk < 0: return
    rng = np.random.default_rng(seed + pi)
    S = rng.integers(0, P, size=(T, len(rows))).astype(np.float32)
    print(f"h={h} D={D} pi={pi} chunk {chunk}/{nchunks}: {len(mine)} of {len(mus)} monomials, rows {len(rows)}, T={T}; "
          f"b-ops {time.time()-t0:.0f}s", flush=True)
    cols = {}; bad = 0; leak = 0; m_rows = n ** 3
    mask = np.ones(m_rows, dtype=bool); mask[rows] = False
    for im, mu in enumerate(mine):
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
        Fs = [Fe(e, fp, -2 * h) if e <= 0 else Fe(e, fm, -wm) for e in es]
        if not np.array_equal(Fe(0, fp, -2 * h), Fe(0, fm, -wm)): bad += 1
        A = np.stack(Fs, axis=1)
        leak += int(A[mask].any())
        Pm = proj(S, A[rows], P)
        for i, e in enumerate(es): cols[(mu, e)] = Pm[:, i].copy()
        if im % 3 == 2:
            from appb_exact import apply_mono; apply_mono.cache_clear()
        if im % 10 == 0: print(f"  {im}/{len(mine)} {time.time()-t0:.0f}s", flush=True)
    out = dict(cols=cols, bad=bad, leak=leak)
    if chunk == 0:
        rhs = []
        for k0 in range(n):
            if lp[k0] != pi: continue
            v = np.zeros((n, n, n), dtype=np.int64); v[k0] = np.eye(n, dtype=np.int64)
            v = v.reshape(-1); assert not v[mask].any()
            rhs.append(v[rows])
        rp = proj(S, np.stack(rhs, axis=1), P); out['rhs'] = [rp[:, i].copy() for i in range(len(rhs))]
    pickle.dump(out, open(f"part_h{h}_D{D}_p{pi}_T{T}_{chunk}of{nchunks}.pkl", "wb"))
    print(f"done chunk {chunk}: two-route failures {bad}, parity leaks {leak}, {time.time()-t0:.0f}s", flush=True)

if __name__ == "__main__":
    if sys.argv[1] == "merge":
        h, D, pi, T, nch = [int(a) for a in sys.argv[2:7]]
        cols = {}; rhs = None; bad = leak = 0
        for c in range(nch):
            d = pickle.load(open(f"part_h{h}_D{D}_p{pi}_T{T}_{c}of{nch}.pkl", "rb"))
            cols.update(d['cols']); bad += d['bad']; leak += d['leak']
            if 'rhs' in d: rhs = d['rhs']
        from compress import compress
        data = dict(h=h, D=D, c=None, P=P, cols=cols, rhs=rhs, n=None, T=T)
        out = compress(data); out['n'] = 0; out['c'] = 0
        print(f"merged h={h} D={D} pi={pi}: {len(cols)} cols, failures {bad}, leaks {leak}, r_all={out['r_all']} (T={T})")
        pickle.dump(out, open(f"Fpar_h{h}_D{D}_p{pi}_T{T}_cmp.pkl", "wb"))
    else:
        h, D, pi, T, chunk, nch = [int(a) for a in sys.argv[1:7]]
        part(h, D, pi, T, chunk, nch)
