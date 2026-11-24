#!/usr/bin/env python3
"""G14e-1: single crossing T_N[x,y,u,v]=[x=v][y=u].  Two-term trade-off
   T = A + B,  rk_row(A)=rank of A as (x,y)|(u,v) <= r,  rk_col(B)=rank of B as (x,u)|(y,v) <= c.
   Equivalent: SWAP_N = X + Y, operator-Schmidt-rank(X) <= r, matrix rank(Y) <= c.
   Numerical feasibility by Levenberg-Marquardt from random starts (real field)."""
import sys, itertools, numpy as np
from scipy.optimize import least_squares

def crossing(N):
    T = np.zeros((N,)*4)
    for x in range(N):
        for y in range(N):
            T[x, y, y, x] = 1
    return T

def feasible(N, r, c, tries=20, seed=0, tol=1e-9, verbose=False):
    T = crossing(N); n = N*N
    rng = np.random.default_rng(seed)
    def unpack(p):
        i = 0
        U = p[i:i+n*r].reshape(n, r); i += n*r
        V = p[i:i+n*r].reshape(n, r); i += n*r
        X = p[i:i+n*c].reshape(n, c); i += n*c
        Y = p[i:i+n*c].reshape(n, c)
        return U, V, X, Y
    def res(p):
        U, V, X, Y = unpack(p)
        A = (U @ V.T).reshape(N, N, N, N)                       # (x,y),(u,v)
        B = (X @ Y.T).reshape(N, N, N, N).transpose(0, 2, 1, 3)  # (x,u),(y,v) -> x,y,u,v
        return (A + B - T).ravel()
    best = 1e9; bestp = None
    for t in range(tries):
        p0 = rng.normal(size=2*n*(r+c)) * 0.7
        s = least_squares(res, p0, method='trf' if 2*n*(r+c) > n*n else 'lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
        f = np.sum(s.fun**2)
        if f < best: best, bestp = f, s.x
        if f < tol: break
    return best, bestp

if __name__ == "__main__":
    N = int(sys.argv[1]); tries = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    n = N*N
    print(f"N={N}: floor r+c>=N={N}; trivial r+c=N^2={n}; eps-type bound 1+N(N-1)/2={1+N*(N-1)//2}")
    for r in range(0, n+1):
        # smallest c feasible
        for c in range(0, n+1):
            if r + c < N: continue
            if r+c > 1+N*(N-1)//2 + 1: break
            f, _ = feasible(N, r, c, tries=tries, seed=100*r+c)
            print(f"  r={r} c={c}: best residual^2 = {f:.3e}", "FEASIBLE" if f < 1e-9 else "", flush=True)
            if f < 1e-9: break
        if r > 1+N*(N-1)//2: break
