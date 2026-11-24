"""G18d: big-M MILP on the DUAL (quadrature) form: min |supp w| s.t. sum_e w_e ev_e = lambda on W_{D'}, |w_e|<=U.
Upper bounds are verified exactly (dual rank mod p + primal consistency mod p). The MILP dual bound is a lower
bound ONLY under |w_e| <= U (after column scaling) and for THIS tau: not a proof.
usage: dmilp.py l d kind mode tl [U]"""
import sys, time, numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
import scipy.sparse as sp
from quad import *
l, d, kind, mode, tl = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4], float(sys.argv[5])
U = float(sys.argv[6]) if len(sys.argv) > 6 else 30.0
tau = taus_for(l, kind); M, bs, nd = mat_float(l, d, tau); ne, nb = M.shape
sc = np.abs(M).max(axis=0); Ms = M / sc; lam = np.zeros(nb); lam[0] = 1
nv = 2 * ne + (1 if mode == "max" else 0)
rows = [sp.hstack([sp.csr_matrix(Ms.T), sp.csr_matrix((nb, nv - ne))])]; lo = list(lam); hi = list(lam)
I = sp.identity(ne); Z1 = sp.csr_matrix((ne, nv - 2 * ne))
rows.append(sp.hstack([I, -U * I, Z1])); lo += [-np.inf] * ne; hi += [0] * ne
rows.append(sp.hstack([-I, -U * I, Z1])); lo += [-np.inf] * ne; hi += [0] * ne
cost = np.zeros(nv); cost[ne:2 * ne] = 1
if mode == "max":
    for k in range(l):
        r = np.zeros(nv); r[[ne + e for e in range(ne) if nd[e][1] == k]] = 1; r[-1] = -1
        rows.append(sp.csr_matrix(r)); lo.append(-np.inf); hi.append(0)
    cost[ne:2 * ne] = 1e-3; cost[-1] = 1
A = sp.vstack(rows).tocsr()
lb = np.concatenate([-U * np.ones(ne), np.zeros(nv - ne)]); ub = np.concatenate([U * np.ones(ne), np.ones(ne), [ne] * (nv - 2 * ne)])
integ = np.concatenate([np.zeros(ne), np.ones(nv - ne)])
t0 = time.time()
res = milp(cost, constraints=LinearConstraint(A, lo, hi), integrality=integ, bounds=Bounds(lb, ub), options={"time_limit": tl})
if res.x is None: print(f"({l},{d}) {kind} {mode}: no solution status {res.status}"); sys.exit()
B = [nd[e] for e in range(ne) if res.x[ne + e] > .5]
rk, ok = exact_dual(l, d, tau, B); pr = exact_primal(l, d, tau, B) if l <= 8 else None
print(f"(l,d)=({l},{d}) tau={kind} mode={mode} U={U}: |B|={len(B)} per_dir={per_dir(l,B)} max={max(per_dir(l,B))} milp_status={res.status} "
      f"dual_bound={res.mip_dual_bound:.3f} exact_dual_ok={ok} rank={rk} primal_ok={pr} bound={bound(l,d)} layered={len(layered(l,d))} "
      f"max|w|={np.abs(res.x[:ne]).max():.2f} {time.time()-t0:.0f}s", flush=True)
print("   B =", sorted(B))
