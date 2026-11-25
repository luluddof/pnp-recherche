# G25-P4: all modes e^g_(n) between weight spaces of the coset module M = V_{alpha/2+L} (A1^n, xi = alpha/2, shifted grading)
# mode(g, Ein, Eout) = matrix [out][in] of the unique mode of e^g mapping M_Ein -> M_Eout  (n = Ein - Eout - 1 since wt e^g = 0)
import sys
sys.path.insert(0, '/home/user/pnp-recherche/code/gen24/g24f')
from objs import Copy, Fr, rankQ
from trank import O_apply, wt
import numpy as np

class Setup:
    def __init__(self, n, E0):
        self.n = n; self.E0 = E0
        G = [[2 if i == j else 0 for j in range(n)] for i in range(n)]
        self.G = G
        self.C = Copy(G, (Fr(1, 2),)*n, coset=[Fr(1, 2)]*n)
        self.S = self.C.S
        self.Mb = {E: self.C.basis(Fr(E) - Fr(n, 4), True) for E in range(E0 + 1)}
        self.idx = {E: {b: i for i, b in enumerate(self.Mb[E])} for E in self.Mb}
        self.A = self.C.basis(0)  # V'_0 basis (e^g, g in S)
        self._m = {}
    def mode(self, g, Ein, Eout):
        key = (tuple(g), Ein, Eout)
        if key in self._m: return self._m[key]
        Min = self.Mb[Ein]; Mout = self.Mb[Eout]; idx = self.idx[Eout]; r = self.n
        M = [[Fr(0)]*len(Min) for _ in Mout]
        eo = Fr(Eout) - Fr(self.n, 4)
        for j, (p, mono) in enumerate(Min):
            q = tuple(p[c] + g[c] for c in range(r)); k2 = eo - wt(self.G, self.C.xi, q)
            if k2 < 0 or k2.denominator != 1: continue
            for m2, c2 in O_apply(self.G, tuple(g), p, mono, int(k2)).items():
                M[idx[(q, m2)]][j] += c2
        self._m[key] = M
        return M

def npf(M): return np.array([[float(x) for x in r] for r in M]) if M and M[0] else np.zeros((len(M), 0))
