import sys, time, json
from shape import *
N, q, p, c = map(int, sys.argv[1:5])
t0 = time.time(); cells = sat_shape(N, q, p, c)
print(f"N={N} (p,q)=({p+c},{q+c}) [S={q},T={p},U={c}]:", "SAT" if cells else "UNSAT", f"{time.time()-t0:.1f}s", flush=True)
if cells:
    show(N, cells); ok, pp, qq = verify_shape(N, cells)
    json.dump({"N": N, "p": pp, "q": qq, "cells": cells}, open(f"shape_{N}_{pp}_{qq}_S{q}T{p}U{c}.json", "w"))
