#!/usr/bin/env python3
"""G19a c2: the 'simultaneous correlation' problem behind A'' for r = 2.
Wire types c in {-1,0,1}^u mod +-, weights w_c (fraction of wires).  A''/M = V(w) := max_{x_c in [-1,1]} min_i |sum_c w_c x_c c_i|
 = max_{s in +-1^u} min_{lambda in simplex} sum_c w_c |<lambda*s, c>|.   m_i = sum_c w_c |c_i| (fraction cut by i).
Designer minimises V(w)/min_i m_i.  Question: is the minimum ~ 1/u (plank bound) or ~ 1/sqrt(u)?
mode 'pm': types in {+-1}^u (every wire cut by every direction, m_i = 1);  mode 'tern': types in {-1,0,1}^u \ 0."""
import itertools, sys, random, time
import numpy as np
from scipy.optimize import linprog

def types(u, mode):
    out = []
    for c in itertools.product((1, 0, -1) if mode == 'tern' else (1, -1), repeat=u):
        nz = [x for x in c if x]
        if nz and nz[0] == 1: out.append(c)
    return np.array(out, dtype=float)

def V(w, C, u):
    best = -1
    for s in itertools.product((1, -1), repeat=u - 1):
        s = np.array((1,) + s, dtype=float)
        # max t : sum_c w_c x_c s_i c_i >= t ; variables x (k), t
        k = len(w); A = np.hstack([-(C * s).T * w, np.ones((u, 1))]); b = np.zeros(u)
        res = linprog(np.r_[np.zeros(k), -1.0], A_ub=A, b_ub=b, bounds=[(-1, 1)] * k + [(None, None)], method='highs')
        best = max(best, -res.fun)
    return best

def score(w, C, u):
    m = (np.abs(C) * w[:, None]).sum(0).min()
    return V(w, C, u) / m if m > 1e-9 else 9.0

if __name__ == '__main__':
    mode = sys.argv[1]; u = int(sys.argv[2]); budget = float(sys.argv[3]); rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 0)
    C = types(u, mode); k = len(C); t0 = time.time()
    wu = np.ones(k) / k
    print(f"mode={mode} u={u} types={k}: uniform weights V/m = {score(wu, C, u):.4f}   (1/u={1/u:.4f}, 0.564/sqrt u={0.564/u**.5:.4f}, 1/sqrt u={u**-.5:.4f})", flush=True)
    best = 9; bw = None
    while time.time() - t0 < budget:
        # random sparse start then local search
        supp = rng.choice(k, size=rng.integers(2, min(k, 2 * u) + 1), replace=False)
        w = np.zeros(k); w[supp] = rng.random(len(supp)) + 0.2; w /= w.sum(); cur = score(w, C, u)
        step = 0.3; stale = 0
        while step > 1e-3 and time.time() - t0 < budget:
            w2 = w.copy(); j = rng.integers(k)
            w2[j] = max(0.0, w2[j] + step * rng.normal() * (w.max()))
            if w2.sum() <= 0: continue
            w2 /= w2.sum(); sc = score(w2, C, u)
            if sc < cur - 1e-9: w, cur, stale = w2, sc, 0
            else:
                stale += 1
                if stale > 40: step /= 2; stale = 0
        if cur < best - 1e-6:
            best, bw = cur, w.copy()
            print(f"  [{time.time()-t0:5.0f}s] best V/min m = {best:.4f}  (u*best={u*best:.3f}, sqrt(u)*best={u**.5*best:.3f})  support: "
                  + "; ".join(f"{tuple(int(x) for x in C[j])}:{bw[j]:.3f}" for j in np.argsort(-bw) if bw[j] > 1e-3), flush=True)
    print(f"RESULT mode={mode} u={u}: min found V/min_i m_i = {best:.4f}; u*V = {u*best:.3f}; sqrt(u)*V = {u**.5*best:.3f}")
