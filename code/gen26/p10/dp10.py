# G26-P10: independent DP of the Lemma-U level model (same model as P2/P5, own code): solvers (x=log2 K, e per N).
# level: x' = xa+lta+xb+ltb, e' = max_th min(ea(cos th + 2^-lh + za), eb(sin th + 2^-lg + zb)) + 2^-min(xa,xb)
# bases: plain sweep s(K)/K^2 (x0<=16) and row-sweep (R3-type ordering) solvers on g x h grids.
import math, sys, time, json
import numpy as np
import zlib10 as Z
base = sys.argv[1]; zname = sys.argv[2]; ratio = float(sys.argv[3]) if len(sys.argv) > 3 else 1.02
XMAX = 2.0**15
def zfun(lt, x):
    if zname == 'F': return Z.zF(lt, x)
    if zname == 'R2': return Z.zR2(lt, x)
    if zname == 'S2half': return Z.zS2(lt, x, 0.5)
    if zname == 'S2': return Z.zS2(lt, x, 1.0)
    raise SystemExit('zeta?')
def lev(ea, eb, al, be):
    # f(th) = min(ea(cos th + al), eb(sin th + be)); first term decreasing, second increasing on [0, pi/2].
    # crossing: ea cos th - eb sin th = eb be - ea al  <=>  R cos(th + phi) = eb be - ea al, R = |(ea, eb)|, phi = atan2(eb, ea)
    ea, eb, al, be = np.broadcast_arrays(*(np.asarray(v, float) for v in (ea, eb, al, be)))
    R = np.hypot(ea, eb); phi = np.arctan2(eb, ea)
    th = np.arccos(np.clip((eb*be - ea*al)/R, -1, 1)) - phi
    th = np.clip(th, 0, math.pi/2)
    v = ea*(np.cos(th) + al)                      # = second term at the crossing
    v = np.minimum(v, eb*(np.sin(th) + be))
    v = np.where(ea*(1+al) <= eb*be, ea*(1+al), v)   # first term below second already at th=0
    v = np.where(ea*al >= eb*(1+be), eb*(1+be), v)   # first term above second even at th=pi/2
    return v
def sweepcost(xr, lg, cA, Gam):     # cut <= U g/(4(g-1)) + cA N^2/U_pad + Gam N + 2/g, padding U_pad >= U chosen optimally
    g = 2.0**lg; v = max(xr, 2*math.sqrt(cA))
    return xr*g/(4*(g-1)) + (v - xr)/4 + cA/v + Gam + 2/g
def rowsolver(lg, lh, T):
    cg, Gg = T(lg); ch, Gh = T(lh); g, h = 2.0**lg, 2.0**lh
    ths = np.linspace(0, math.pi/2, 20001)
    def sc(xr, l, cA, Gam):
        gg = 2.0**l; v = np.maximum(xr, 2*math.sqrt(cA))
        return xr*gg/(4*(gg-1)) + (v - xr)/4 + cA/v + Gam + 2/gg
    a = sc(1/h + np.cos(ths), lg, cg, Gg); b = sc(1/g + np.sin(ths), lh, ch, Gh)
    return float(np.minimum(a, b).max()) + 1e-7     # grid max + small margin (Lipschitz, 20001 points)
cache = {}
def T(lg):
    if lg not in cache: cache[lg] = Z.r3terms(lg, base)
    return cache[lg]
cand = []
def ddc(K): return 0 if K <= 1 else (K//2)*(K - K//2) + ddc(K - K//2)
for x0 in range(1, 17):
    K = 2**x0; cand.append((float(x0), (ddc(K) if base == 'poly' else ((K+1)**2)//4 - 1)/K**2))
for lg in (np.arange(6, 40.01, 0.25) if base != 'poly' else []):
    for d in (0.0, 0.5, 1.0, 2.0, 4.0):
        cand.append((2*lg + d, rowsolver(lg, lg + d, T)))
nb = int(math.log(XMAX/0.5)/math.log(ratio)) + 2
bx = np.full(nb, np.nan); be_ = np.full(nb, np.nan); bmu = np.zeros(nb)
def put(xs, es):
    xs = np.ravel(xs); es = np.ravel(es); ok = (xs < XMAX) & (es > 0) & np.isfinite(es)
    for x, e in zip(xs[ok], es[ok]) if xs[ok].size < 50 else []:
        pass
    xs, es = xs[ok], es[ok]
    if xs.size == 0: return
    mu = 1/(es*es*xs); b = np.clip((np.log(xs/0.5)/math.log(ratio)).astype(int), 0, nb-1)
    order = np.argsort(-mu, kind='stable'); b, mu, xs, es = b[order], mu[order], xs[order], es[order]
    _, first = np.unique(b, return_index=True)
    b, mu, xs, es = b[first], mu[first], xs[first], es[first]
    m = mu > bmu[b]; bmu[b[m]] = mu[m]; bx[b[m]] = xs[m]; be_[b[m]] = es[m]
put([c[0] for c in cand], [c[1] for c in cand])
lts = np.arange(1.0, 48.01, 1.0)
zc = {}
def zrow(x):
    if x not in zc: zc[x] = np.array([zfun(lt, x) for lt in lts])
    return zc[x]
t0 = time.time()
for i in range(nb):
    if np.isnan(bx[i]): continue
    xa, ea = bx[i], be_[i]; za = zrow(xa)
    js = [j for j in range(i + 1) if not np.isnan(bx[j]) and xa + bx[j] + 2 < XMAX]
    if not js: continue
    xb = bx[js][:, None, None]; eb = be_[js][:, None, None]
    zb = np.array([zrow(bx[j]) for j in js])[:, None, :]
    lg = xa + lts[None, :, None]; lh = xb + lts[None, None, :]
    e2 = lev(ea, eb, 2.0**-lh + za[None, :, None], 2.0**-lg + zb) + 2.0**-np.minimum(xa, xb)
    put(np.broadcast_to(lg + lh, e2.shape), e2)
ok = ~np.isnan(bx); X = bx[ok]; E = be_[ok]; o = np.argsort(X); X, E = X[o], E[o]; Em = np.minimum.accumulate(E)
json.dump({'x': X.tolist(), 'e': E.tolist()}, open(f'sol10_{base}_{zname}_{ratio}.json', 'w'))
print(f"base={base} zeta={zname} ratio={ratio} time {time.time()-t0:.0f}s")
for a in (2, 4):
    s = []
    for Y in (1000, 2000, 4000, 8000):
        ys = np.exp(np.linspace(math.log(Y), math.log(2*Y), 600)); k = np.searchsorted(X, ys, side='right') - 1
        s.append(f"[{Y},{2*Y}] {(2*Em[k]*np.sqrt(a*ys)).max():.3f}")
    print(f"  a={a}: " + "; ".join(s))
