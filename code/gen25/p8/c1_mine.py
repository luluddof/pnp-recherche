import math, sys
sys.path.insert(0, '.')
from recur_mine import *
# c1 = limsup 2 e_l sqrt(a x_{l+1}); schedules: zeta target z0*q^l, also z0*q^l with other q; x0 grid.
for name, base, split, entry, a in [("C' brute,split,entry1,a=9", 'brute', True, 1, 9),
                                    ("C  brute,paper-zeta,entry2,a=9", 'brute', False, 2, 9),
                                    ("D  sweep,a=9", 'sweep', True, 1, 9),
                                    ("E1 sweep,a=4", 'sweep', True, 1, 4),
                                    ("E2 sweep,a=2", 'sweep', True, 1, 2)]:
    best = None
    for x0 in range(4, 16):
        for k in range(2, 10):
            for q in (0.45, 0.5, 0.55, 0.6, 0.7):
                v, xs, es = c1_limsup(base, split, entry, a, x0, lambda l: 2.0**-k * q**l, L=34)
                lim = max(v[24:]); 
                if best is None or lim < best[0]: best = (lim, x0, k, q, v, xs, es)
    lim, x0, k, q, v, xs, es = best
    print(f"{name:32s} c1~{lim:.3f} x0={x0} z0=2^-{k} q={q}  v[l=4,8,16,24,33]={[round(v[i],3) for i in (4,8,16,24,33)]}  liminf ratio={2*es[33]*math.sqrt(a*xs[33])/v[33]:.4f}", flush=True)
