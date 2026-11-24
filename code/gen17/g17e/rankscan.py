import sys
from quad import *
for a in sys.argv[1:]:
    l, d = map(int, a.split(",")); tau = taus_for(l, "rand"); L = layered(l, d); rk, ok = exact_dual(l, d, tau, L)
    print(f"(l,d)=({l},{d}) D'={l-d-1}: bound={bound(l,d)} layered={len(L)} rank(pruned)={rk} feasible={ok} pruned/bound={rk/bound(l,d):.2f} layered/bound={len(L)/bound(l,d):.2f}", flush=True)
