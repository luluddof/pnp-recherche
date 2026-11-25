# Row-sweep base: check P13's q (geometric delta grid 1e-5..1e3, 160 pts, chord bound via monotonicity) against
# chords on a much finer grid incl. (0, 1e-5]: delta in [1e-10, 1e3], 3000 pts. Own max-min via bisection (sc increasing).
import numpy as np, math, zlib10 as Z
def sc(xr, l, cA, Gam):
    gg = 2.0**l; v = np.maximum(xr, 2*math.sqrt(cA))
    return xr*gg/(4*(gg-1)) + (v - xr)/4 + cA/v + Gam + 2/gg
def G(d, lg, lh, cg, Gg, ch, Gh):
    g, h = 2.0**lg, 2.0**lh; dp = np.minimum(d, 1.0)
    fa = lambda th: sc(((1-dp)*np.cos(th) + dp)/(1+d) + 1/h, lg, cg, Gg)
    fb = lambda th: sc(((1-dp)*np.sin(th) + dp)/(1+d) + 1/g, lh, ch, Gh)
    lo = np.zeros_like(d); hi = np.full_like(d, math.pi/2)
    for _ in range(100):
        mid = (lo + hi)/2; up = fa(mid) > fb(mid); lo = np.where(up, mid, lo); hi = np.where(up, hi, mid)
    return (1 + d)*np.maximum(np.minimum(fa(lo), fb(lo)), np.minimum(fa(hi), fb(hi)))
DQ = np.concatenate([[0.0], np.geomspace(1e-5, 1e3, 160)])
DF = np.geomspace(1e-10, 1e3, 3000)
worst = 0
for lg in (6, 8, 10, 12, 16, 20, 24, 30):
    for dd in (0.0, 1.0, 4.0):
        lh = lg + dd; cg, Gg = Z.r3terms(lg, 'R3pp'); ch, Gh = Z.r3terms(lh, 'R3pp')
        Gq = G(DQ, lg, lh, cg, Gg, ch, Gh); G0 = Gq[0]
        q13 = max(float(((Gq[2:] - G0 + 1e-12)/(DQ[1:-1]*G0)).max()), float((Gq[1] - G0 + 1e-12)/(DQ[1]*G0)))
        Gf = G(DF, lg, lh, cg, Gg, ch, Gh); qf = float(((Gf - G0)/(DF*G0)).max()); ksmall = float(((Gf[DF < 1e-5] - G0)/(DF[DF < 1e-5]*G0)).max())
        worst = max(worst, qf - q13)
        print(f"lg={lg} dd={dd}: q13={q13:.4f}  fine-grid sup chord={qf:.4f}  sup on (1e-10,1e-5)={ksmall:.4f}", flush=True)
print(f"max(fine - q13) = {worst:.2e}")
