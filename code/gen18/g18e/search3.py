import numpy as np
from search import solve
wires = [(0, 7), (1, 6), (2, 5), (3, 4)]; u = 3
for prof in [(3,1,1),(2,2,1)]:
    for seed in range(4):
        for nf in (250, 1000, 3000):
            try:
                c, x, t = solve(wires, u, prof, tries=1, seed=seed, maxnfev=nf)
                print(f"pk3m4 {prof} seed {seed} max_nfev {nf}: residual^2 = {c:.3e}  max|param| = {np.abs(x).max():.3e}", flush=True)
            except Exception as e:
                print(f"pk3m4 {prof} seed {seed} max_nfev {nf}: solver error {type(e).__name__}", flush=True)
