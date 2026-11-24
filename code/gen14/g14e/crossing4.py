import sys, numpy as np
from crossing import feasible
N = int(sys.argv[1]); pairs = eval(sys.argv[2]); tries = int(sys.argv[3])
for (r, c) in pairs:
    f, _ = feasible(N, r, c, tries=tries, seed=7*r+c)
    print(f"N={N} r={r} c={c}: best residual^2 = {f:.3e}", "FEASIBLE" if f < 1e-9 else "not found", flush=True)
