"""Random search over binary [D,k] codes containing all-ones (k=4,5): best achievable
error-rank fraction (translation-invariant UB, exact MILP) at each gain threshold, vs the cube."""
import sys, numpy as np
sys.path.insert(0, '.')
from crossing_design import code_structure, ti_rank_binary, words
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
k = int(sys.argv[1]); K = 2 ** k
thresholds = [0.2, 0.25, 1/3, 0.4, 0.5]
best = {g: (1.0, None) for g in thresholds}
def evalG(G, tag):
    phi, C, w = code_structure(G)
    D = G.shape[1]
    for d in range(D):
        r, _ = ti_rank_binary(G, w, d)
        for g in thresholds:
            if d / D <= g + 1e-9 and r / K < best[g][0]:
                best[g] = (r / K, tag)
        if r == 1: break
evalG(np.eye(k, dtype=np.int64), f"cube Q{k}")
for D in range(k + 1, 3 * k + 1):
    for trial in range(25 if k == 4 else 8):
        G = rng.integers(0, 2, size=(k, D)); G[0, :] = 1
        Z = words(k)
        if len({tuple(r) for r in (Z @ G % 2)}) < K: continue
        evalG(G, f"[{D},{k}] #{trial}")
for g in thresholds:
    print(f"k={k} gain<={g:.3f}: best rank fraction {best[g][0]:.4f} by {best[g][1]}")
