#!/usr/bin/env python3
"""G18f: hexagon (2x3 grid) two-term problem in operator form.
P = S_01 (x) S_02 (x) S_12 on 3 sites, each site = 2 qudits of dim N (slot toward each other site).
B = triangle operator: sum_{al,be,ga} B0[o0,i0,al,ga] B1[o1,i1,al,be] B2[o2,i2,be,ga]; A = P - B, want rank(A) <= r.
Objective f_r(B) = sum_{i>r} sigma_i(P - B)^2  (A eliminated exactly by Eckart-Young)."""
import numpy as np, itertools
from scipy.optimize import minimize

def perm_P(N):
    d = N*N
    # site i slots: site0=(->1,->2), site1=(->0,->2), site2=(->0,->1)
    P = np.zeros((d,)*6)
    for a, a2, b, b2, c, c2 in itertools.product(range(N), repeat=6):
        # in: site0=(a,a2) site1=(b,b2) site2=(c,c2);  swaps: a<->b (01), a2<->c (02), b2<->c2 (12)
        i0, i1, i2 = a*N+a2, b*N+b2, c*N+c2
        o0, o1, o2 = b*N+c, a*N+c2, a2*N+b2
        P[o0, o1, o2, i0, i1, i2] = 1.0
    return P.reshape(d**3, d**3)

def build(B0, B1, B2):
    d = B0.shape[0]
    M = np.einsum('aiAG,bjAB,ckBG->abcijk', B0, B1, B2, optimize=True)
    return M.reshape(d**3, d**3)

def tail(R, r):
    U, s, Vt = np.linalg.svd(R, full_matrices=False)
    f = float((s[r:]**2).sum())
    G = (U[:, r:] * s[r:]) @ Vt[r:]          # R - R_r
    return f, G, s

def make_obj(P, d, D, r):
    D01, D12, D02 = D
    sh = [(d, d, D01, D02), (d, d, D01, D12), (d, d, D12, D02)]
    sz = [int(np.prod(s)) for s in sh]
    def unpack(x):
        o = 0; out = []
        for s, n in zip(sh, sz):
            out.append(x[o:o+n].reshape(s)); o += n
        return out
    def fg(x):
        B0, B1, B2 = unpack(x)
        R = P - build(B0, B1, B2)
        f, G, s = tail(R, r)
        G6 = (-2*G).reshape((d,)*6)
        g0 = np.einsum('abcijk,bjAB,ckBG->aiAG', G6, B1, B2, optimize=True)
        g1 = np.einsum('abcijk,aiAG,ckBG->bjAB', G6, B0, B2, optimize=True)
        g2 = np.einsum('abcijk,aiAG,bjAB->ckBG', G6, B0, B1, optimize=True)
        return f, np.concatenate([g0.ravel(), g1.ravel(), g2.ravel()])
    return fg, unpack, sum(sz)

def search(P, d, D, r, starts=10, maxiter=600, seed=0, scale=0.5, verbose=False):
    fg, unpack, n = make_obj(P, d, D, r)
    rng = np.random.default_rng(seed)
    best = (np.inf, None); vals = []
    for s in range(starts):
        x0 = rng.normal(size=n) * scale
        res = minimize(fg, x0, jac=True, method='L-BFGS-B', options=dict(maxiter=maxiter, maxfun=maxiter*2, ftol=1e-15, gtol=1e-10))
        vals.append(res.fun)
        if res.fun < best[0]: best = (res.fun, res.x.copy())
        if res.fun < 1e-9: break
    return best[0], best[1], sorted(vals)
