# dim span{ a[-h]|_{V_h} : wt a <= w }  (functionals on V_h coming from lower-weight states)
import sys; sys.dont_write_bytecode = True
import numpy as np, flint
from fcomp import Ops
from appb_exact import partitions
P = 16777213
for h in range(4, int(sys.argv[1]) + 1):
    O = Ops(P); n = len(O.basis(h)); out = []
    for w in range(0, h + 1):
        rows = [O.mat(mu, -h, h).reshape(-1) for ww in range(w + 1) for mu in partitions(ww)]
        M = flint.nmod_mat(len(rows), n, (np.array(rows) % P).ravel().tolist(), P)
        out.append(M.rank())
    print(f"h={h} p(h)={n}: dim by max weight w=0..h:", out, flush=True)
