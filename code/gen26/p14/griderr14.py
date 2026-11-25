# base e: P10/P13 use a 20001-point theta grid max + 1e-7; compare with the exact max-min (bisection). Under-estimate?
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
THS = np.linspace(0, math.pi/2, 20001); worst = 0
for lg in np.arange(6, 40.01, 0.25):
    for dd in (0.0, 0.5, 1.0, 2.0, 4.0):
        lh = lg + dd; cg, Gg = Z.r3terms(lg, 'R3pp'); ch, Gh = Z.r3terms(lh, 'R3pp'); g, h = 2.0**lg, 2.0**lh
        e = float(np.minimum(sc(1/h + np.cos(THS), lg, cg, Gg), sc(1/g + np.sin(THS), lh, ch, Gh)).max()) + 1e-7
        G0 = float(G(np.array([0.0]), lg, lh, cg, Gg, ch, Gh)[0]); worst = max(worst, (G0 - e)/G0)
print(f"max relative under-estimate of base e by grid (after +1e-7 margin): {worst:.2e}")
