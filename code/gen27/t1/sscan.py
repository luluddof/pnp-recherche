# G27-T1: k_max(h, D) via mode-span reduction (direct operator form, P1 direct_cols), parity sectors, exact rows, mod P.
# S_r(D)  = span{ t^mu[r] : V_h -> V_{h+r},  |mu| <= D }   (rho side, r in [-h, Wmax-h])
# S'_r(D) = span{ t^mu[r] : V_{h-r} -> V_h,  |mu| <= D }   (kappa side, r in [h-Wpmax, h])
# columns: rho: nu -> t^nu[-r] A (A in S_r, r <= W-h);  kappa: nu -> B t^nu[-r] (B in S'_r, r >= h-W').
# feasible(W,W') iff  e_x^* (x) Id in span for all x.  Parity (alpha -> -alpha) splits everything into 2 sectors.
# D = 'inf': S_r = full parity-graded Hom (density of modes).
# usage: python3 sscan.py h Dlist Wmax Wpmax [P] [Wlist]
import sys, time, pickle, json
sys.dont_write_bytecode = True
import numpy as np, flint
from appb_exact import partitions, apply_mono
from fcomp import Ops, mm

def fl(A, P): return flint.nmod_mat(A.shape[0], A.shape[1], (A % P).ravel().tolist(), P)
def rowbasis(A, P):
    if A.shape[0] == 0: return A
    R, r = fl(A, P).rref()
    return np.array(R.entries(), dtype=np.int64).reshape(A.shape)[:r]
def rank(A, P):
    if A.shape[0] == 0: return 0
    return fl(A, P).rank()

class Setup:
    def __init__(self, h, Wmax, Wpmax, P):
        self.h, self.P = h, P; self.O = O = Ops(P); self.Vh = Vh = O.basis(h); n = self.n = len(Vh)
        self.Wmax, self.Wpmax = Wmax, Wpmax
        self.rr = list(range(-h, Wmax - h + 1)); self.rk = [r for r in range(h - Wpmax, h + 1)]
        par = lambda E: np.array([len(s) % 2 for s in O.basis(E)], dtype=np.int64)
        self.par = par
        lp = par(h)
        self.rowpar = ((lp[:, None, None] + lp[None, :, None] + lp[None, None, :]) % 2).reshape(-1)
        self.rows = {pi: np.nonzero(self.rowpar == pi)[0] for pi in (0, 1)}
        # b-operators
        self.Bm = {r: np.stack([O.mat(nu, -r, h + r) for nu in Vh]) for r in self.rr}     # (nu, sigma, tau) tau in V_{h+r}
        self.Bp = {r: np.stack([O.mat(nu, -r, h) for nu in Vh]) for r in self.rk}         # (nu, tau, sigma') tau in V_{h-r}
        # parity masks for Hom spaces
        self.hmask = {}
        for r in self.rr:
            pt = par(h + r); self.hmask[('rho', r)] = ((pt[:, None] + lp[None, :]) % 2).reshape(-1)   # [tau, sigma']
        for r in self.rk:
            pt = par(h - r); self.hmask[('kap', r)] = ((lp[:, None] + pt[None, :]) % 2).reshape(-1)   # [sigma, tau]
    def cols(self, side, r, basis, pi):
        """basis: (k, q*n) or (k, n*q') -> column vectors (k, n^3) restricted to sector rows"""
        h, n, P = self.h, self.n, self.P; k = basis.shape[0]
        if k == 0: return np.zeros((0, len(self.rows[pi])), dtype=np.int64)
        if side == 'rho':
            q = self.Bm[r].shape[2]
            As = basis.reshape(k, q, n).transpose(1, 0, 2).reshape(q, k * n)
            C = mm(self.Bm[r].reshape(n * n, q), As, P).reshape(n, n, k, n).transpose(2, 0, 1, 3)
        else:
            q = self.Bp[r].shape[1]
            Bs = basis.reshape(k * n, q)
            Bpa = self.Bp[r].transpose(1, 0, 2).reshape(q, n * n)
            C = mm(Bs, Bpa, P).reshape(k, n, n, n).transpose(0, 2, 1, 3)
        return C.reshape(k, n ** 3)[:, self.rows[pi]]
    def rhs(self, pi):
        n = self.n; lp = self.par(self.h); out = []
        for x in range(n):
            if lp[x] != pi: continue
            v = np.zeros((n, n, n), dtype=np.int64); v[x] = np.eye(n, dtype=np.int64); out.append(v.reshape(-1)[self.rows[pi]])
        return np.array(out, dtype=np.int64)

def scan(st, S, label, Ws, log):
    """S[(side, r, pi)] = basis rows. Returns {W: (minW', k)}"""
    h, P = st.h, st.P; t0 = time.time()
    res = {}
    for pi in (0, 1):
        R = st.rhs(pi)
        cr = {r: rowbasis(st.cols('rho', r, S[('rho', r, pi)], pi), P) for r in st.rr}
        ck = {r: rowbasis(st.cols('kap', r, S[('kap', r, pi)], pi), P) for r in st.rk}
        # prefix spans
        Rho = {}; acc = np.zeros((0, len(st.rows[pi])), dtype=np.int64)
        for W in range(0, st.Wmax + 1):
            r = W - h
            if r in cr: acc = rowbasis(np.concatenate([acc, cr[r]]), P)
            Rho[W] = acc
        Kap = {-1: np.zeros((0, len(st.rows[pi])), dtype=np.int64)}; acc = Kap[-1]
        for Wp in range(0, st.Wpmax + 1):
            r = h - Wp
            if r in ck: acc = rowbasis(np.concatenate([acc, ck[r]]), P)
            Kap[Wp] = acc
        def feas(W, Wp):
            M = np.concatenate([Rho[W], Kap[Wp]]); a = rank(M, P); b = rank(np.concatenate([M, R]), P); return a == b, b - a
        for W in Ws:
            if W > st.Wmax: continue
            hi = min(4 * h - 1 - W, st.Wpmax)
            if not feas(W, hi)[0]: res.setdefault(W, []).append(None); continue
            lo = -1
            while lo < hi:
                mid = (lo + hi) // 2
                if feas(W, mid)[0]: hi = mid
                else: lo = mid + 1
            res.setdefault(W, []).append(lo)
        print(f"   [{label}] sector {pi} done, rank Rho(Wmax)={Rho[st.Wmax].shape[0]} Kap(Wpmax)={Kap[st.Wpmax].shape[0]} "
              f"rows {len(st.rows[pi])}  {time.time()-t0:.0f}s", flush=True)
    out = {}
    for W, l in res.items():
        if None in l: out[W] = None
        else: Wp = max(l); out[W] = (Wp, 4 * h - 1 - W - Wp)
    line = f"h={h} D={label}: " + " ".join(f"W{W}:{'-' if out[W] is None else str(out[W][0])+'/k'+str(out[W][1])}" for W in sorted(out))
    print(line, flush=True); log.write(line + "\n"); log.flush()
    return out

import os
ADJ = os.environ.get('ADJ', '1') == '1'
def zl(lam):
    import math
    z = 1
    for k in set(lam): m = lam.count(k); z *= k ** m * math.factorial(m)
    return z
def adjoint_fill(st, S, P):
    """S'_r(D) = { A^dagger : A in S_{-r}(D) },  A^dagger = G_h^{-1} A^T G_{h-r}, G = diag(z_lambda) (Shapovalov form;
    Y(a,z)^dagger = Y(e^{zL1}(-z^-2)^{L0} a, 1/z) keeps the weight filtration)."""
    h, n = st.h, st.n
    gh_inv = np.array([pow(zl(s), P - 2, P) for s in st.O.basis(h)], dtype=np.int64)
    for r in st.rk:
        src = st.O.basis(h - r); q = len(src); g = np.array([zl(s) % P for s in src], dtype=np.int64)
        for pi in (0, 1):
            A = S[('rho', -r, pi)]                       # (k, q*n) [tau, sigma'] : V_h -> V_{h-r}
            k = A.shape[0]
            B = A.reshape(k, q, n).transpose(0, 2, 1)    # A^T : (k, n, q) [sigma, tau]
            B = (B * gh_inv[None, :, None]) % P
            B = (B * g[None, None, :]) % P
            S[('kap', r, pi)] = B.reshape(k, n * q)
def main(h, Dlist, Wmax, Wpmax, P=16777213, Ws=None):
    t0 = time.time(); st = Setup(h, Wmax, Wpmax, P); n = st.n
    print(f"h={h} n={n} Wmax={Wmax} W'max={Wpmax} P={P}: setup {time.time()-t0:.0f}s", flush=True)
    Ws = Ws or list(range(0, Wmax + 1))
    keys = [('rho', r, pi) for r in st.rr for pi in (0, 1)] + [('kap', r, pi) for r in st.rk for pi in (0, 1)]
    def fulldim(key): return int((st.hmask[key[:2]] == key[2]).sum())
    S = {}
    for key in keys:
        q = st.O.basis(h + key[1]) if key[0] == 'rho' else st.O.basis(h - key[1])
        S[key] = np.zeros((0, len(q) * n), dtype=np.int64)
    log = open(f"sscan_h{h}_P{P}.txt", "a"); results = {}
    finite = [D for D in Dlist if D != 'inf']
    import os
    ck = f"Sck_h{h}_P{P}_W{Wmax}_{Wpmax}.pkl"; w0 = 0
    if os.path.exists(ck):
        S, w0 = pickle.load(open(ck, 'rb')); w0 += 1; print(f" resumed from checkpoint, weights < {w0} done", flush=True)
    for w in range(w0, (max(finite) if finite else -1) + 1):
        new = {key: [] for key in keys}
        for mu in partitions(w):
            pi = len(mu) % 2
            for r in st.rr:
                if S[('rho', r, pi)].shape[0] < fulldim(('rho', r, pi)): new[('rho', r, pi)].append(st.O.mat(mu, r, h).ravel())
            if not ADJ:
              for r in st.rk:
                if S[('kap', r, pi)].shape[0] < fulldim(('kap', r, pi)): new[('kap', r, pi)].append(st.O.mat(mu, r, h - r).ravel())
        for key in keys:
            if new[key]: S[key] = rowbasis(np.concatenate([S[key], np.array(new[key], dtype=np.int64) % P]), P)
        if ADJ: adjoint_fill(st, S, P)
        pickle.dump((S, w), open(ck + '.tmp', 'wb')); os.replace(ck + '.tmp', ck)
        dims = sum(S[k].shape[0] for k in keys); full = sum(fulldim(k) for k in keys)
        nfull = sum(S[k].shape[0] == fulldim(k) for k in keys)
        msg = f" weight {w}: total mode-span dim {dims} / {full} (full Hom); saturated blocks {nfull}/{len(keys)}  {time.time()-t0:.0f}s"
        print(msg, flush=True); log.write(msg + "\n")
        if w in finite:
            results[w] = scan(st, S, str(w), Ws, log)
        if dims == full:
            print(" mode span saturated: all larger D equal D=inf", flush=True); log.write(f" saturated at weight {w}\n")
            for D in finite:
                if D > w: results[D] = results.get(w) or scan(st, S, str(w), Ws, log)
            break
    if 'inf' in Dlist:
        Sf = {}
        for key in keys:
            m = st.hmask[key[:2]]; idx = np.nonzero(m == key[2])[0]
            E = np.zeros((len(idx), len(m)), dtype=np.int64); E[np.arange(len(idx)), idx] = 1; Sf[key] = E
        results['inf'] = scan(st, Sf, 'inf', Ws, log)
    pickle.dump(results, open(f"sscan_h{h}_P{P}_W{Wmax}_{Wpmax}.pkl", "wb"))
    print(f"total {time.time()-t0:.0f}s")

if __name__ == "__main__":
    h = int(sys.argv[1]); Dlist = [a if a == 'inf' else int(a) for a in sys.argv[2].split(',')]
    Wmax, Wpmax = int(sys.argv[3]), int(sys.argv[4]); P = int(sys.argv[5]) if len(sys.argv) > 5 else 16777213
    Ws = [int(a) for a in sys.argv[6].split(',')] if len(sys.argv) > 6 else None
    main(h, Dlist, Wmax, Wpmax, P, Ws)
