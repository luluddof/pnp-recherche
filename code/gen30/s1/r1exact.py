# R(k,1) exactly: min over Paulis P=X^a Z^b, real kappa of rank(D_f - kappa P) = D - max real-eigenvalue geometric multiplicity of P^{-1} D_f.
# Twisted translations map (a,b) -> (a, b+Ja), J invertible, so b = 0 WLOG. P^{-1} D_f = X^a D_f (X^a involution): a signed permutation.
import numpy as np, sys
from core import *
for k in (1, 2, 3, 4):
    D = 4**k; f = beta(k); z = np.arange(D); best = (D, None)
    # X^a D_f |z> = f(z)|z+a>: 2-cycles {z, z+a} (a != 0) with block [[0, f(z+a)],[f(z), 0]] -> eigenvalues +-sqrt(f(z)f(z+a))
    for a in range(D):
        if a == 0:
            mult = max((f == 1).sum(), (f == -1).sum())
        else:
            prod = f * f[z ^ a]               # per 2-cycle (counted twice)
            npos = (prod == 1).sum() // 2      # cycles with eigenvalues +1,-1 (real)
            mult = npos                        # eigenvalue +1 multiplicity (also -1); prod=-1 cycles give +-i (not real)
        r = D - mult
        if r < best[0]: best = (r, a)
    Q = lambda a: popc((a & ((1 << k) - 1)) & (a >> k)) & 1
    print(f'k={k} D={D}: R(1) = {best[0]} (a={best[1]}), Bell p(1) = {D - max((f==1).sum(), (f==-1).sum())}, line D/4 = {D/4}', flush=True)
