import numpy as np, sys
from wreath2 import setup
from wreath import mat
for ns in [[2,2,2],[3,3],[4,4],[2,2,2,2]]:
    S, SW, E, w, LR, sR, LC, rng = setup(ns)
    r = S['r']
    out = []
    for f, wt in zip(E, w):
        Fs = np.array([(f @ mat(np.eye(r)[o], S) @ f).ravel() for o in range(r)])
        m2 = np.linalg.matrix_rank(Fs, tol=1e-7)
        m = int(round(np.sqrt(m2))); sg = int(round(np.trace(SW @ f) / np.trace(f)))
        out.append((wt, m, wt // max(m,1), sg))
    print(ns, "pieces (weight, mult m, dim d, swap sign):", out, flush=True)
