# Interpolation in one torus-weight space: degree-d polys with leg contents (cx,cy,cu,cv) vanishing on J(p,q) (random points mod PR),
# and whether some of them do not vanish on Z3 = {M(x,v) M'(y,u)} (the closure of the swap's orbit).
# Probabilistic-exact: ranks mod PR of evaluation matrices at (#monomials + 40) random points.
import numpy as np, sys, time, json
from flint import nmod_mat
from common import *
from monos import monos, parts
def evalmat(pts, mon):
    P = np.array([X.ravel() for X in pts], dtype=np.int64)
    E = np.ones((len(pts), len(mon)), dtype=np.int64)
    for k in range(mon.shape[1]): E = (E * P[:, mon[:, k]]) % PR
    return E
def nrank(E):
    return nmod_mat(E.shape[0], E.shape[1], E.ravel().tolist(), PR).rank()
def test(cont, N=3, p=2, q=2, seed=0, extra=40):
    mon = np.array(monos(cont), dtype=np.int64); m = len(mon)
    rng = np.random.default_rng(seed)
    EJ = evalmat([ptJ(N,p,q,rng) for _ in range(m + extra)], mon)
    rJ = nrank(EJ)
    if rJ == m: return m, 0, 0
    EZ = evalmat([ptZ3(N,rng) for _ in range(40)], mon)
    rJZ = nrank(np.vstack([EJ, EZ]))
    return m, m - rJ, rJZ - rJ          # (#monomials, dim of equations of J in this weight, dim of their image on Z3)
if __name__ == "__main__":
    dmax = int(sys.argv[1]); mmax = int(sys.argv[2]); mode = sys.argv[3] if len(sys.argv) > 3 else 'sym'
    p, q = (int(sys.argv[4]), int(sys.argv[5])) if len(sys.argv) > 5 else (2, 2)
    for d in range(1, dmax+1):
        P = parts(d)
        if mode == 'sym': W = [(mu,nu,nu,mu) for mu in P for nu in P if nu <= mu]
        else:  # all dominant weights up to the D4 x leg symmetries (crudely: all, dedup by sorting is not valid in general)
            W = [(a,b,c,e) for a in P for b in P for c in P for e in P]
        for cont in W:
            t0 = time.time()
            n = len(monos(cont))
            if n > mmax: print(d, cont, n, 'SKIPPED', flush=True); continue
            m, a, b = test(cont, p=p, q=q)
            print(d, cont, f"monomials={m} eqs_of_J={a} nonzero_on_Z3={b}", 'WITNESS' if b > 0 else '', f"{time.time()-t0:.0f}s", flush=True)
