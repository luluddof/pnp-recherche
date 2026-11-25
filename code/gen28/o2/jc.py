# Jordan-Chevalley: M_row = R S^-1 = M_s + M_n. Is (R_s = M_s S, C_s = S - R_s) again a decomposition, and what is osr(C_s)?
import sys, sympy as sp, numpy as np
from pieces import *
for f in sys.argv[1:]:
    N, (R, C, s) = load_json(f); S = s*swap(N)
    M = sp.Matrix(rowflat(R).tolist()) * sp.Matrix(rowflat(S).tolist()).inv()
    P, J = M.jordan_form()
    D = sp.diag(*[J[i, i] for i in range(J.shape[0])])
    Ms = P * D * P.inv()
    Rs = Ms * sp.Matrix(rowflat(S).tolist())
    n = N*N
    Rs_t = np.array(Rs.tolist(), dtype=object).reshape(N, N, N, N)
    Cs_t = (S.astype(object) - Rs_t)
    Ccol = sp.Matrix(Cs_t.transpose(0, 2, 1, 3).reshape(n, n).tolist())
    print(f.split('/')[-1], "p =", sp.Matrix(rowflat(R).tolist()).rank(), " q =", sp.Matrix(colflat(C).tolist()).rank(),
          "| rank R_s =", Rs.rank(), " osr C_s =", Ccol.rank(), " entries rational:", all(x.is_rational for x in Rs), flush=True)
