# for each orbit rep support S (3 Paulis), numeric min over kappa in R^3 of rank(D_f - sum kappa P): AP on fixed support,
# plus a sigma-minimisation cross-check (Nelder-Mead on the r+1-th singular value) for the target r = 5.
import numpy as np, pickle, sys, time
from core import *
from sl import trunc
from scipy.optimize import minimize
k = 2; D = 16; V = np.diag(beta(k)); reps = pickle.load(open('k2q3_reps.pkl', 'rb')); pts = pickle.load(open('k2_pts.pkl', 'rb'))
rng = np.random.default_rng(0); T0 = time.time()
def mats(S): return [pauli(pts[i][0], pts[i][1], 4) for i in S]
def feasible_AP(Ps, r, starts=8):
    M = np.array([P.ravel() for P in Ps]).T; Q, _ = np.linalg.qr(M)
    for s in range(starts):
        L = trunc(rng.standard_normal((D, D)) * 2, r)
        for it in range(3000):
            R = (V - L).ravel(); C = (Q @ (Q.T @ R)).reshape(D, D); L = trunc(V - C, r)
            if it % 200 == 199 and np.linalg.norm(V - C - L) < 1e-10: return True
        if np.linalg.norm(V - C - L) < 1e-8: return True
    return False
def sigmin(Ps, r, starts=30):
    best = 9
    f = lambda x: np.linalg.svd(V - sum(xi * P for xi, P in zip(x, Ps)), compute_uv=False)[r]
    for s in range(starts):
        x0 = rng.standard_normal(len(Ps)) * 2
        res = minimize(f, x0, method='Nelder-Mead', options={'xatol': 1e-12, 'fatol': 1e-14, 'maxiter': 4000})
        best = min(best, res.fun)
    return best
hist = {}; worst = []
for i, S in enumerate(reps):
    Ps = mats(S); rmin = None
    for r in range(4, 9):
        if feasible_AP(Ps, r): rmin = r; break
    hist[rmin] = hist.get(rmin, 0) + 1
    if rmin is not None and rmin <= 5: worst.append((S, rmin))
    if rmin is None or rmin >= 6:
        sm = sigmin(Ps, 5, 12)
        if sm < 1e-6: worst.append((S, 'sigma6~0', sm))
print('rank histogram over orbit reps (None = >8):', hist, f'{time.time()-T0:.0f}s')
print('reps with rank <= 5 or sigma_6 -> 0:', worst)
