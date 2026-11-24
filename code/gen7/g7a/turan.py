"""Compare: (a) combinatorial min edge set hitting all Q_{d+1} subcubes (total and per-direction max) via MILP;
(b) whether such optimal hitting sets are algebraically realizable."""
import sys, numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from core import *
def hit(l, m, mode="total", tl=120):
    E = edges(l); C = subcubes(l, m); ne = len(E)
    nv = ne + (1 if mode == "max" else 0)
    A = np.zeros((len(C) + (l if mode == "max" else 0), nv))
    for i, es in enumerate(C): A[i, es] = 1
    lo = [1] * len(C); hi = [np.inf] * len(C)
    c = np.zeros(nv)
    if mode == "max":
        for k in range(l):
            r = len(C) + k
            for e, (z, w, kk) in enumerate(E):
                if kk == k: A[r, e] = 1
            A[r, -1] = -1; lo.append(-np.inf); hi.append(0)
        c[-1] = 1
    else: c[:] = 1
    integ = np.ones(nv)
    ub = np.ones(nv); 
    if mode == "max": ub[-1] = np.inf
    res = milp(c, constraints=LinearConstraint(A, lo, hi), integrality=integ, bounds=Bounds(0, ub),
               options={"time_limit": tl})
    bad = [e for e in range(ne) if res.x[e] > .5]
    return res, bad, E
if __name__ == "__main__":
    l, d = int(sys.argv[1]), int(sys.argv[2]); mode = sys.argv[3] if len(sys.argv) > 3 else "total"
    res, bad, E = hit(l, d + 1, mode)
    per = [sum(1 for e in bad if E[e][2] == k) for k in range(l)]
    good = [E[e] for e in range(len(E)) if e not in set(bad)]
    print(f"l={l} d={d} {mode}: comb-hitting optimum total={len(bad)} per={per} status={res.status} "
          f"dual={res.mip_dual_bound:.2f} alg-feasible(generic tau)={feasible(l,d,good)}")
