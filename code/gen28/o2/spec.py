# Spectral structure of M_row = R_row S_row^{-1} and M_col = C_col S_col^{-1} (S scaled so R + C = S).
# If M is diagonalizable, the k-th power type-class piece A_i (row side) is diagonal in the eigen-word basis with entry
#   coefficient of T^i in prod_s (m_s + T (1 - m_s)),  so its rank is a count of words (no rank subadditivity loss beyond that).
import sys, numpy as np, sympy as sp
from pieces import *
def analyse(name, N, R, C, s):
    S = s*swap(N)
    out = {}
    for side, X, Sf in (('row', rowflat(R), rowflat(S)), ('col', colflat(C), colflat(S))):
        M = sp.Matrix(X.tolist()) * sp.Matrix(Sf.tolist()).inv()
        ev = M.eigenvals()            # exact
        diag = all((M - l*sp.eye(M.shape[0])).rank() == (M - l*sp.eye(M.shape[0]))**2 and True or True for l in [])
        mins = {}
        for l, m in ev.items():
            A = M - l*sp.eye(M.shape[0]); r1 = A.rank(); r2 = (A*A).rank()
            mins[l] = (m, M.shape[0] - r1, r1 == r2)     # (alg mult, geom mult, semisimple)
        out[side] = mins
        print(f"{name} {side}: eigen (value: alg mult, geom mult, semisimple) = {mins}", flush=True)
    return out
if __name__ == '__main__':
    for which in sys.argv[1:]:
        if which == 'oct': N, (R, C, s) = octonion()
        elif which.startswith('F1_'): N, (R, C, s) = F1(int(which[3:]))
        else: N, (R, C, s) = load_json(which)
        analyse(which, N, R, C, s)
