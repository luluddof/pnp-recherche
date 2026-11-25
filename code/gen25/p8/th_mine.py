import math, sys
sys.path.insert(0, '.')
from recur_mine import *
rows = [("C  paper-zeta, entry 2/K, brute, paper top", 'paper9', 'brute', False, 2),
        ("C' split, entry 1/K, brute, paper top", 'paper9', 'brute', True, 1),
        ("D  + sweep", 'paper9', 'sweep', True, 1),
        ("E1 + T4 top (det.)", 'E1', 'sweep', True, 1),
        ("E2 + 8' top (LV)", 'E2', 'sweep', True, 1)]
for r in [float(a) for a in sys.argv[1:]] or [4.26, 1.0]:
    for tau, tn in [(1.0, '2^n'), (math.log2(1.6423), '1.6423^n'), (math.log2(1.30704), '1.30704^n')]:
        for name, var, base, split, entry in rows:
            res = best_threshold(r, tau, var, base=base, split=split, entry=entry)
            th, x0, lts, xs, es = res
            print(f"r={r} {tn:10s} {name:42s} log2 m* = {th:9.1f}  x0={x0} levels={len(lts)} lt={[round(a,2) for a in lts]} log2K_top={xs[-1]:.1f} e_top={es[-1]:.5f}", flush=True)
