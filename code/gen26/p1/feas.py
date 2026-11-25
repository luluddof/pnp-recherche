# Feasibility of the degree-1 ansatz for given (W, W'):  l (x) Id in span{F^mu_e : e allowed}, for all l.
# allowed (mu,e):  e <= W' - 2h   (kappa side, d_r)   or   e >= |mu| + h - W   (rho side, c_r).
import sys, time, pickle
sys.dont_write_bytecode = True
import numpy as np, flint

def allowed(h, mu, e, W, Wp):
    return e <= Wp - 2 * h or e >= sum(mu) + h - W

def rank_mod(cols, P):
    if not cols: return 0
    A = np.stack(cols, axis=1) % P
    M = flint.nmod_mat(A.shape[0], A.shape[1], A.ravel().tolist(), P)
    return M.rank()

def feasible(data, W, Wp, Dmax=None, extra=None):
    h, P = data['h'], data['P']
    cols = [v for (mu, e), v in data['cols'].items() if allowed(h, mu, e, W, Wp) and (Dmax is None or sum(mu) <= Dmax)]
    ra = rank_mod(cols, P); rb = rank_mod(cols + data['rhs'], P)
    return ra == rb, ra, rb, len(cols)

def direct_cols(h, D, P, W, Wp, G):
    """validation: columns computed directly, O^-_{mu,r} r in [-h, W-h], O^+_{mu,r} r in [h-W', h]"""
    from fcomp import Ops, mm
    from appb_exact import partitions
    O = Ops(P); Vh = O.basis(h); n = len(Vh); c = G.shape[1]
    mus = [m for w in range(D + 1) for m in sorted(partitions(w))]
    cols = []
    for mu in mus:
        for r in range(-h, W - h + 1):
            Au = mm(O.mat(mu, r, h), G, P)
            v = np.stack([mm(O.mat(nu, -r, h + r), Au, P) for nu in Vh])      # [nu, sigma, j]
            cols.append(v.reshape(-1))
        for r in range(h - Wp, h + 1):
            v = np.stack([mm(O.mat(mu, r, h - r), mm(O.mat(nu, -r, h), G, P), P) for nu in Vh])
            cols.append(v.reshape(-1))
    return cols

if __name__ == "__main__":
    fn = sys.argv[1]; mode = sys.argv[2] if len(sys.argv) > 2 else "scan"
    data = pickle.load(open(fn, "rb")); h = data['h']
    print(f"h={h} D={data['D']} c={data['c']} rows={data['n']**2*data['c']} cols={len(data['cols'])}", flush=True)
    if mode == "validate":
        for W in range(0, 3 * h + 1):
            for Wp in range(0, 3 * h + 1):
                if W + Wp > 4 * h or W + Wp < 2 * h: continue
                f1 = feasible(data, W, Wp)
                dc = direct_cols(h, data['D'], data['P'], W, Wp, data['G'])
                ra = rank_mod(dc, data['P']); rb = rank_mod(dc + data['rhs'], data['P'])
                print(f"W={W} W'={Wp} k={4*h-1-W-Wp}: F-route feasible={f1[0]} (ranks {f1[1]},{f1[2]}) | direct feasible={ra==rb} ({ra},{rb})", flush=True)
    else:
        # for each W, minimal W' (binary search; monotone in W')
        t0 = time.time()
        for W in range(0, 4 * h):
            lo, hi = 0, 4 * h - 1 - W   # threshold W' = 4h-1-W is feasible (App B)
            if hi < 0: break
            while lo < hi:
                mid = (lo + hi) // 2
                if feasible(data, W, mid)[0]: hi = mid
                else: lo = mid + 1
            f = feasible(data, W, lo)
            print(f"W={W:>3}: min W'={lo:>3}  k = {4*h-1-W-lo}   (feasible check {f[0]}, rank {f[1]}, cols {f[3]})  {time.time()-t0:.0f}s", flush=True)

def probe(data, Ws, ks, Dmax=None):
    h = data['h']; out = []
    for W in Ws:
        for k in ks:
            Wp = 4 * h - 1 - W - k
            if Wp < 0: continue
            t = time.time(); f = feasible(data, W, Wp, Dmax)
            print(f"  W={W:>3} W'={Wp:>3} k={k}: feasible={f[0]} ranks ({f[1]},{f[2]}) cols {f[3]}  [{time.time()-t:.0f}s]", flush=True)
            out.append((W, Wp, k, f[0], f[1], f[2]))
    return out
