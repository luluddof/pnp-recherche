# E2: required truncation t for a correct decision on random 3-SAT, Delta=4.26
import random, time, math, sys
from math import comb
import numpy as np
from ie import *
D = float(sys.argv[1]) if len(sys.argv) > 1 else 4.26
rng = random.Random(2026)
print(f"Delta={D}")
print(" n   m  inst  #UNSAT  mean_rmax  rmax/n  rmax/m  mean_tdec(hist)  tdec/n  log2C(m,<=t)/n")
rows = []
for n in [10, 14, 18, 22, 26, 30, 35, 40]:
    m = round(D * n); R = []; TD = []; U = 0; L = []
    K = 12
    for i in range(K):
        F = rand3sat(n, m, rng)
        sat = is_sat(n, F); U += (not sat)
        rm = rmax_maxsat(n, F); R.append(rm)
        if n <= 22:
            a = hist(n, F); assert max(i for i in range(len(a)) if a[i]) == rm
            td, ans = t_decide(a); assert (ans == 'SAT') == sat; TD.append(td)
            if not sat: assert td == rm + (rm % 2)       # theorem: UNSAT needs smallest even t >= rmax
        else:
            # UNSAT: t_dec = rmax rounded up to even (proved). SAT: use rmax as proxy (upper bound+1)
            TD.append(rm + (rm % 2) if not sat else rm)
        t = TD[-1]
        L.append(math.log2(sum(comb(m, k) for k in range(t + 1))) / n)
    mr = np.mean(R); mt = np.mean(TD)
    rows.append((n, mr, mt))
    print(f"{n:3d} {m:3d} {K:4d} {U:6d}   {mr:8.2f}  {mr/n:6.3f}  {mr/m:6.3f}   {mt:8.2f}{'' if n<=22 else '*'}       {mt/n:6.3f}   {np.mean(L):6.2f}")
    sys.stdout.flush()
ns = np.array([r[0] for r in rows]); ts = np.array([r[2] for r in rows])
A = np.vstack([ns, np.ones_like(ns)]).T
c, d = np.linalg.lstsq(A, ts, rcond=None)[0]
p = np.polyfit(np.log(ns), np.log(ts), 1)
print(f"fit t ~ {c:.3f} n + {d:.2f};  power fit t ~ n^{p[0]:.3f}")
print("(* = n>22: tdec from exact MaxSAT rmax; exact for UNSAT, proxy for SAT)")
