# Lower bound for ALL two-piece blockwise schemes at block size lambda0 (eigen-type model):
# enumerate every spectral class tuple (z,k,m) for row and column, subject only to the
# rank constraint P(N) <= rank_col(R)+rank_row(C) and <= rank_row(R)+rank_col(C) with P(N) >= N.
import numpy as np, math
from rho import f_curve, TAUS
def fast_f(cls, taus=TAUS, G=201):
    z, k, m = cls
    lw = lambda w: math.log2(w) if w > 0 else None
    A = np.linspace(0, 1, G)
    best = np.full(len(taus), -np.inf)
    for a in A:
        if a > 0 and z == 0: continue
        r = 1 - a
        base = (-a*math.log2(a) if a > 0 else 0) + (a*math.log2(z) if a > 0 else 0)
        for i, t in enumerate(taus):
            if a > t + 1e-12: continue
            lo = max(t - a, 0.0)
            if lo > r + 1e-12: continue
            if m == 0 and lo > 1e-12: continue
            if k == 0: pm = r
            elif m == 0: pm = 0.0
            else: pm = max(lo, r*m/(m+k))
            pm = min(pm, r); pk = r - pm
            v = base; bad = False
            for p, w in ((pm, m), (pk, k)):
                if p > 1e-15:
                    if w == 0: bad = True; break
                    v += -p*math.log2(p) + p*math.log2(w)
            if bad: continue
            best[i] = max(best[i], v)
    return best
for lam, PN in ((1, 2), (2, 4), (3, 8)):
    N2 = 4**lam
    tuples = [(z, k, N2-z-k) for z in range(N2+1) for k in range(N2+1-z)]
    taus = np.linspace(0, 1, 101)
    F = {c: fast_f(c, taus) for c in tuples}
    bestv, arg = 9, None
    for rc in tuples:
        for cc in tuples:
            zr, kr, mr = rc; zc, kc, mc = cc
            if (N2 - zc) + (N2 - kr) < PN or (N2 - zr) + (N2 - kc) < PN: continue
            if min(N2 - zr, N2 - zc, N2 - kr, N2 - kc) < 1: continue  # R, C nonzero
            v = np.max(F[rc] + F[cc]) / (4*lam)
            if v < bestv: bestv, arg = v, (rc, cc)
    print(f"lambda0={lam}: model lower bound for 2-piece blockwise rho >= {bestv:.4f}  attained by row{arg[0]} col{arg[1]}")
print("asymptotic (lambda0->inf, exponent LP): rho >= (2 + c*)/4 >= 3/4, c* = lim log P(N)/log N")
