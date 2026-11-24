import sys, itertools, numpy as np
from search import solve
wires = [(0, 7), (1, 6), (2, 5), (3, 4)]; u = 3
for prof in [(2,2,0),(2,1,1),(4,0,0),(3,1,0),(3,1,1),(2,2,1),(3,2,0),(2,2,2),(3,2,1),(3,3,0),(4,1,1),(4,2,0)]:
    best = np.inf
    for seed in range(6):
        try:
            c, x, t = solve(wires, u, prof, tries=1, seed=seed, maxnfev=250)
            best = min(best, c)
        except Exception as e:
            pass
        if best < 1e-18: break
    print(f"pk3m4 profile {prof} N={sum(prof)}: best residual^2 = {best:.3e}", "FOUND" if best < 1e-18 else "", flush=True)
