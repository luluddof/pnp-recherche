# Asymptotic constant c1 in  time <= 2^{c1 W / sqrt(log2 m)}  for the re-optimised recursion.
# e_l from the exact recursion (thmA_explicit.levels, rho_bar -> 1/sqrt2 as lam -> inf), x_l = log2 K_l.
# Level l is usable once log2 m >= a x_l + b; the worst m for level l is just below the start of level l+1, so
# c1 = limsup_l 2 e_l sqrt(a x_{l+1} + b)   (top-level deviation -> 0 inside the range, as m/K^a -> inf).
import math, thmA_explicit as T
def c1_of(x0, z0, q, base, split, entry1, a, b, L=26):
    ts = T.sched_ts(x0, L, z0, q, 'hoeff'); lv = T.levels(x0, ts, base, split, entry1, T.RHO)
    vals = [2*lv[l][1]*math.sqrt(a*lv[l+1][0] + b) for l in range(L)]
    return vals
for name, base, split, e1, a, b in [("paper recursion as written, a=9 (W<=3m)", None, None, None, 9, 0),
                                   ("C' exact Hoeffding+split, brute base, a=9", 'brute', True, True, 9, 0),
                                   ("D  + sweep base, a=9", 'sweep', True, True, 9, 0),
                                   ("E1 + Chebyshev top, a=4", 'sweep', True, True, 4, 0),
                                   ("E2 + Las Vegas top, a=2", 'sweep', True, True, 2, 0)]:
    if base is None:
        e, x = 0.5, 18; vals = []; xs=[x]; es=[e]
        for l in range(26):
            z = 2.0**(-l-5); e = (T.RHO+3*z)*e + (3/2**x if x < 1000 else 0); x = 2*(x + 3*l + 18); xs.append(x); es.append(e)
        vals = [2*es[l]*math.sqrt(9*xs[l+1]) for l in range(26)]
        print(f"{name:45s} 2e_l sqrt(a x_(l+1)), l=10,15,20,25: {[round(vals[l],2) for l in (10,15,20,25)]}  (x W<=3m: per m multiply by 3)")
        continue
    best = None
    for x0 in range(4, 16):
        for k in range(2, 10):
            for q in (0.5, 0.6, 0.7, 0.8):
                v = c1_of(x0, 2.0**-k, q, base, split, e1, a, b)
                lim = max(v[18:])
                if best is None or lim < best[0]: best = (lim, x0, 2.0**-k, q, v)
    lim, x0, z0, q, v = best
    print(f"{name:45s} c1 = limsup = {lim:.3f}  (x0={x0}, z0={z0}, q={q});  l=4,6,8,10,15,20,25: {[round(v[l],3) for l in (4,6,8,10,15,20,25)]}")
