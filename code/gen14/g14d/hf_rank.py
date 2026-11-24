# G14d: how many INDEPENDENT conditions does a layered bad set impose on the dual space W_{D'} (D' = l-d-1)?
# Commutative class, one point per direction. W_{D'} = span{ chi_S(z) pi_S(t) t^i : |S|+i <= D' }, ev_(z,k)(W) = W(z,tau_k).
# Exact duality: B feasible  <=>  lambda (W -> sum_z W(z,0)) lies in span{ev_e : e in B}.  Counting bounds can reach at most rank(ev_B).
import itertools, numpy as np
from math import comb
p = 2147483629
def rank_mod(M):
    M = M.copy() % p; rk = 0; rows, cols = M.shape
    for c in range(cols):
        nz = np.nonzero(M[rk:, c])[0]
        if len(nz) == 0: continue
        i = rk + nz[0]
        if i != rk: M[[rk, i]] = M[[i, rk]]
        inv = pow(int(M[rk, c]), p-2, p)
        M[rk] = (M[rk]*inv) % p
        nzr = np.nonzero(M[:, c])[0]
        for r in nzr:
            if r != rk: M[r] = (M[r] - M[r, c]*M[rk]) % p
        rk += 1
        if rk == rows: break
    return rk
def run(l, d):
    Dp = l-d-1
    if Dp < 0: return
    taus = list(range(1, l+1))
    inv = [pow(t, p-2, p) for t in taus]
    basis = [(S, i) for s in range(Dp+1) for S in itertools.combinations(range(l), s) for i in range(Dp-s+1)]
    best = None
    for shift in range(d+1):
        cuts = sorted(set([0] + list(range(shift if shift > 0 else d+1, l+1, d+1)) + [l+1]))
        bnd = cuts[1:-1]                       # boundary between weight b-1 and b
        B = [(z, k) for z in range(1 << l) for k in range(l) if not z >> k & 1 and (bin(z).count('1')+1) in bnd]
        if best is None or len(B) < len(best): best = B
    B = best
    def row(z, k, with_lambda=False):
        r = []
        t = taus[k]
        for S, i in basis:
            if k in S: r.append(0); continue
            v = pow(t, i, p)
            for j in S:
                v = v*((1 - t*inv[j]) % p) % p
                if z >> j & 1: v = (-v) % p
            r.append(v)
        return r
    M = np.array([row(z, k) for z, k in B], dtype=np.int64)
    lam = np.array([[(1 << l) % p if (len(S) == 0 and i == 0) else 0 for S, i in basis]], dtype=np.int64)
    rk = rank_mod(M); rk2 = rank_mod(np.vstack([M, lam]))
    new = sum(comb(l, j)*max(0, l-d-2*j) for j in range(l+1))
    print(f"l={l} d={d}: |B_layered|={len(B)}  rank(ev_B on W_D')={rk}  lambda in span: {rk2 == rk}  dim W_D'={len(basis)}  new bound={new}", flush=True)
for l, d in [(3,1),(4,1),(4,2),(5,1),(5,2),(5,3),(6,2),(6,3),(7,2),(7,3),(7,4),(8,3),(8,4),(9,3),(9,4),(10,4),(10,5)]:
    run(l, d)
