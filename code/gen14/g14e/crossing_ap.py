#!/usr/bin/env python3
"""Crossing two-term trade-off by LM with analytic Jacobian (fast).  T = A + B, rk_row(A)<=r, rk_col(B)<=c."""
import sys, numpy as np
from scipy.optimize import least_squares
def run(N, r, c, tries, seed=0, tol=1e-10):
    n = N*N; rng = np.random.default_rng(seed)
    T = np.zeros((N,)*4)
    for x in range(N):
        for y in range(N): T[x, y, y, x] = 1
    Trow = T.reshape(n, n)
    perm = np.arange(n*n).reshape(N, N, N, N).transpose(0, 2, 1, 3).reshape(-1)   # col-flat index -> row-flat position
    In = np.eye(n)
    def unpack(p):
        i = 0; out = []
        for k in (r, r, c, c):
            out.append(p[i:i+n*k].reshape(n, k)); i += n*k
        return out
    def res(p):
        U, V, X, Y = unpack(p)
        Rm = (U @ V.T).ravel() - Trow.ravel()
        Bc = (X @ Y.T).ravel()
        Rm[perm] += Bc
        return Rm
    def jac(p):
        U, V, X, Y = unpack(p)
        JU = np.kron(In, V)            # d (U V^T)[i,j] / d U[i,k] = V[j,k]
        JV = np.einsum('ik,jl->ijlk', U, In).reshape(n*n, n*r)
        JX = np.zeros((n*n, n*c)); JY = np.zeros((n*n, n*c))
        JX[perm] = np.kron(In, Y); JY[perm] = np.einsum('ik,jl->ijlk', X, In).reshape(n*n, n*c)
        return np.hstack([JU, JV, JX, JY])
    best = 1e9; nfeas = 0; hist = []
    for t in range(tries):
        p0 = rng.normal(size=2*n*(r+c)) * rng.choice([0.3, 0.7, 1.2])
        s = least_squares(res, p0, jac=jac, method='trf', xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=400)
        f = float(np.sum(s.fun**2)); hist.append(round(f, 3))
        if f < best: best, bp = f, s.x
        if f < tol: nfeas += 1; break
    return best, sorted(set(hist))[:4], bp
if __name__ == "__main__":
    N = int(sys.argv[1]); pairs = eval(sys.argv[2]); tries = int(sys.argv[3])
    for (r, c) in pairs:
        best, hist, bp = run(N, r, c, tries, seed=31*r+c)
        print(f"N={N} (r,c)=({r},{c}): best residual^2={best:.2e}  lowest local minima {hist}  ->", "FEASIBLE" if best < 1e-10 else "not found", flush=True)
