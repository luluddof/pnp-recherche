# Compute the Laurent polynomials F^mu_e(nu) (one-variable picture of Lemma 9) for all monomials mu of weight <= D,
# applied to c random vectors gamma (sketch of the sigma' index), all nu (basis of V_h), all sigma.  Arithmetic mod P.
# Output: npz with Fcols[(mu,e)] -> vector over (nu, sigma, j), plus the RHS for l = e_nu0^*, plus checks.
import sys, math, time, pickle
sys.dont_write_bytecode = True
import numpy as np
from math import comb
from appb_exact import partitions, apply_mono

def mm(A, B, P):
    """(A @ B) mod P exactly, A,B int64 with entries in [0,P), P < 2^24, inner dim < 2^13."""
    A = np.ascontiguousarray(A, dtype=np.float64)
    B1 = (B >> 12).astype(np.float64); B0 = (B & 4095).astype(np.float64)
    r1 = np.fmod(A @ B1, P).astype(np.int64); r0 = np.fmod(A @ B0, P).astype(np.int64)
    return ((r1 << 12) + r0) % P

class Ops:
    def __init__(self, P):
        self.P = P; self.bas = {}; self.idx = {}
    def basis(self, E):
        if E not in self.bas:
            b = sorted(partitions(E)) if E >= 0 else []
            self.bas[E] = b; self.idx[E] = {s: i for i, s in enumerate(b)}
        return self.bas[E]
    def mat(self, mu, r, E):
        """matrix of t^mu[r]: V_E -> V_{E+r}, mod P"""
        src = self.basis(E); tgt = self.basis(E + r); ix = self.idx[E + r]
        M = np.zeros((len(tgt), len(src)), dtype=np.int64)
        if not tgt: return M
        for j, s in enumerate(src):
            for s2, c in apply_mono(mu, r, s):
                M[ix[s2], j] = (M[ix[s2], j] + c) % self.P
        return M

def compute(h, D, c, P, seed=1, verbose=True):
    t0 = time.time()
    O = Ops(P)
    Vh = O.basis(h); n = len(Vh)
    mus = [m for w in range(D + 1) for m in sorted(partitions(w))]
    rng = np.random.default_rng(seed)
    G = rng.integers(0, P, size=(n, c), dtype=np.int64)
    # f+ data: O+_{mu}(w) = t^mu[h-w] (V_w->V_h) . t^nu[w-h] (V_h->V_w),  exponent m = w-2h, w in [0,2h]
    # f- data: O-_{mu}(w) = t^nu[h-w] (V_w->V_h) . t^mu[w-h] (V_h->V_w),  exponent m = -w
    BpG = {}; BmS = {}
    wm = max(2 * h, h + D)
    for w in range(0, wm + 1):
        pw = len(O.basis(w))
        if w <= 2 * h: BpG[w] = np.concatenate([mm(O.mat(nu, w - h, h), G, P) for nu in Vh], axis=1)   # (pw, n*c) [nu, j]
        BmS[w] = np.concatenate([O.mat(nu, h - w, w) for nu in Vh], axis=0)            # (n*n, pw) [nu, sigma]
    if verbose: print(f"h={h}: b-operators built {time.time()-t0:.1f}s", flush=True)
    cols = {}; checks = []
    for mu in mus:
        hp = sum(mu); H = h + hp
        fp = {}; fm = {}
        for w in range(0, wm + 1):
          if w <= 2 * h:
            Ad = O.mat(mu, h - w, w)            # V_w -> V_h
            X = mm(Ad, BpG[w], P)               # (n[sigma], n*c [nu,j])
            fp[w - 2 * h] = X.reshape(n, n, c).transpose(1, 0, 2).reshape(-1)   # [nu, sigma, j]
          if True:
            Au = O.mat(mu, w - h, h)            # V_h -> V_w
            Y = mm(BmS[w], mm(Au, G, P), P)     # (n*n [nu,sigma], c)
            fm[-w] = Y.reshape(-1)
        L = n * n * c; Z = np.zeros(L, dtype=np.int64)
        def Fe(e, f, lo, hi):
            acc = np.zeros(L, dtype=np.int64)
            for j in range(0, H + 1):
                m = e - j
                if lo <= m <= hi and m in f:
                    acc = (acc + ((-1) ** j * comb(H, j)) % P * f[m]) % P
            return acc
        Fs = {}
        for e in range(-2 * h, h + hp + 1):
            Fs[e] = Fe(e, fp, -2 * h, 0) if e <= 0 else Fe(e, fm, -wm, 0)
        # check: e = 0 both ways (f- route needs m in [-H, 0] i.e. w <= H <= 2h ok)
        ok0 = np.array_equal(Fe(0, fp, -2 * h, 0), Fe(0, fm, -wm, 0))
        # check: F_e = 0 just outside the support, via the available series
        okout = (not Fe(-2 * h - 1, fp, -2 * h, 0).any())
        if hp + h + 1 <= wm:
            okout = okout and not Fe(h + hp + 1, fm, -2 * h, 0).any()
        checks.append((mu, ok0, okout))
        for e, v in Fs.items(): cols[(mu, e)] = v
    # F(1) check against Gamma: sum_e F_e should be (t^mu[-h] t^nu) * Id
    if verbose:
        print(f"h={h}: F built for {len(mus)} monomials, {len(cols)} columns, {time.time()-t0:.1f}s; "
              f"checks e=0 consistency {all(x[1] for x in checks)}, outside-support {all(x[2] for x in checks)}", flush=True)
    rhs = []
    for k0 in range(n):
        v = np.zeros((n, n, c), dtype=np.int64); v[k0] = G
        rhs.append(v.reshape(-1))
    return dict(h=h, D=D, c=c, P=P, mus=mus, cols=cols, rhs=rhs, G=G, n=n, checks=checks)

if __name__ == "__main__":
    h, D, c = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); P = int(sys.argv[4]) if len(sys.argv) > 4 else 16777213
    res = compute(h, D, c, P)
    # pairing check: sum_e F_e = (t^mu[-h] t^nu) Id
    O = Ops(P); Vh = O.basis(h); n = len(Vh); bad = 0
    for mu in res['mus']:
        tot = sum(res['cols'][(mu, e)] for e in range(-2 * h, h + sum(mu) + 1)) % P
        pair = O.mat(mu, -h, h).reshape(-1) % P   # (1, n): value on t^nu
        exp = np.zeros((n, n, c), dtype=np.int64)
        for k0 in range(n): exp[k0] = pair[k0] * res['G'] % P
        bad += not np.array_equal(tot, exp.reshape(-1))
    print("Lemma 9(c) check sum_e F_e = (t^mu[-h] t^nu) Id: failures", bad)
    pickle.dump(res, open(f"F_h{h}_D{D}_c{c}_P{P}.pkl", "wb"))
