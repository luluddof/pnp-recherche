# ST-only (U=0) SAT with colour-size cap from the bound (S sizes <= g*(p), T sizes <= g*(q)); caps are implied by Thm 3, only help.
import sys, time, json
from shape import *
N, q, p = map(int, sys.argv[1:4])
def gstar(t, N):
    g = 1
    while g + 1 <= N and (g + 1) * g // 2 <= t: g += 1
    return g
ms = max(gstar(p, N), gstar(q, N))
t0 = time.time(); cells = sat_shape(N, q, p, 0, maxsize=ms)
print(f"N={N} (p,q)=({p},{q}) ST-only, maxsize {ms}:", "SAT" if cells else "UNSAT", f"{time.time()-t0:.1f}s", flush=True)
if cells:
    show(N, cells); ok, pp, qq = verify_shape(N, cells)
    json.dump({"N": N, "p": pp, "q": qq, "cells": cells}, open(f"st_{N}_{pp}_{qq}.json", "w"))
