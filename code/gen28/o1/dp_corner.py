# G28-O1: dp13 (G26-P13, mode splitiii) with the App-B circle replaced by a hypothetical frontier
#   F_c = min(circle, c(1-u))  in normalized coordinates (u,v) = (theta/2, T/2)  (one corner at (u,v)=(1,0)),
# or two corners if TWO=1. Levels and row-sweep bases both use the hull H_c; both grid orientations are tried.
# usage: python dp_corner.py c [ratio]   (c = 'inf' -> plain circle, validation against 3.166)
import math, sys, time, json, os
import numpy as np
sys.path.insert(0, '/home/user/pnp-recherche/code/gen26/p13')
import zlib10 as Z
base = 'R3pp'; cs = sys.argv[1]; ratio = float(sys.argv[2]) if len(sys.argv) > 2 else 1.02
TWO = os.environ.get('TWO', '0') == '1'
XMAX = 2.0**float(os.environ.get('LXMAX', '15'))
c = math.inf if cs == 'inf' else float(cs)
if math.isinf(c): uc, vc = 1.0, 0.0
else: uc = (c*c - 1)/(c*c + 1); vc = c*(1 - uc)
if TWO and not math.isinf(c) and c <= 1 + math.sqrt(2): uc = vc = c/(1 + c)   # two segments meet before the circle
phic = math.atan2(vc, uc)
# second corner (mirror) at (0,1): v-side segment from (0,1) to (vc,uc)
phic2 = math.atan2(uc, vc) if TWO else math.pi/2
def hsup(a, b):
    """support function of H_c at (a,b) >= 0."""
    ang = np.arctan2(b, a)
    r = np.hypot(a, b)
    arc = np.where(ang < phic, a*uc + b*vc, np.where(ang > phic2, a*(vc if TWO else 0) + b*(uc if TWO else 1), r))
    return np.maximum(np.maximum(arc, a), b)
def fcov(x, kind):
    if kind == 0: return 1.0
    K = 2.0**min(x, 60)
    if K < 8: return 1.0
    eps = (4*K - 2)/((K - 2)*(K - 3))
    return min(1.0, math.sqrt((1 + eps)/2))
def lev2(ea, eb, al, be):
    """max over the circle of min(ea(u+al), eb(v+be)) and minimax weight w (P13's code)."""
    ea, eb, al, be = np.broadcast_arrays(*(np.asarray(v, float) for v in (ea, eb, al, be)))
    R = np.hypot(ea, eb); phi = np.arctan2(eb, ea)
    th = np.arccos(np.clip((eb*be - ea*al)/R, -1, 1)) - phi
    th = np.clip(th, 0, math.pi/2)
    v = np.minimum(ea*(np.cos(th) + al), eb*(np.sin(th) + be))
    w = eb*np.cos(th)/(ea*np.sin(th) + eb*np.cos(th))
    cA = ea*(1+al) <= eb*be; cB = ea*al >= eb*(1+be)
    v = np.where(cA, ea*(1+al), v); w = np.where(cA, 1.0, w)
    v = np.where(cB, eb*(1+be), v); w = np.where(cB, 0.0, w)
    return v, w
def levH(ea, eb, al, be):
    v, w = lev2(ea, eb, al, be)
    ea, eb, al, be, v, w = np.broadcast_arrays(ea, eb, al, be, v, w); v = v.copy(); w = w.copy()
    if not math.isinf(c):
        DJ = ea*(uc + al) - eb*(vc + be)
        seg = DJ < 0
        us = (eb*(c + be) - ea*al)/(ea + c*eb)
        vs = np.where(us >= 1, ea*(1 + al), ea*(us + al))
        ws = np.where(us >= 1, 1.0, c*eb/(ea + c*eb))
        v = np.where(seg, vs, v); w = np.where(seg, ws, w)
        if TWO:   # mirror corner at (0,1): segment u = c(1-v), v in [uc,1]
            DJ2 = ea*(vc + al) - eb*(uc + be)
            seg2 = DJ2 > 0
            vv = (ea*(c + al) - eb*be)/(eb + c*ea)
            v2 = np.where(vv >= 1, eb*(1 + be), eb*(vv + be))
            w2 = np.where(vv >= 1, 0.0, ea/(ea + c*eb) * 0 + eb/(c*ea + eb))
            v = np.where(seg2, v2, v); w = np.where(seg2, w2, w)
    return v, w
# boundary parametrisation for row-sweep bases: tau in [0, Ltot], tau=0 at (1,0), u decreasing / v increasing
L2s = math.hypot(1 - uc, vc); L1 = phic2 - phic; L3 = math.hypot(vc, 1 - uc) if TWO else 0.0
LT = L2s + L1 + L3
def bpt(tau):
    tau = np.asarray(tau, float)
    s = np.clip(tau/max(L2s, 1e-300), 0, 1)
    us, vs = 1 - s*(1 - uc), s*vc
    ang = phic + np.clip(tau - L2s, 0, L1)
    ua, va = np.cos(ang), np.sin(ang)
    u = np.where(tau <= L2s, us, ua); v = np.where(tau <= L2s, vs, va)
    if TWO:
        s3 = np.clip((tau - L2s - L1)/max(L3, 1e-300), 0, 1)
        u3, v3 = vc*(1 - s3), uc + s3*(1 - uc)
        m = tau > L2s + L1; u = np.where(m, u3, u); v = np.where(m, v3, v)
    return u, v
def sc(xr, l, cA, Gam):
    gg = 2.0**l; v = np.maximum(xr, 2*math.sqrt(cA))
    return xr*gg/(4*(gg-1)) + (v - xr)/4 + cA/v + Gam + 2/gg
DQ = np.concatenate([[0.0], np.geomspace(1e-5, 1e3, 160)])
def rowsolver(lg, lh, T):
    cg, Gg = T(lg); ch, Gh = T(lh); g, h = 2.0**lg, 2.0**lh
    d = DQ; dp = np.minimum(d, 1.0)
    def fa(tau): u, v = bpt(tau); return sc(((1-dp)*u + dp)/(1+d) + 1/h, lg, cg, Gg)
    def fb(tau): u, v = bpt(tau); return sc(((1-dp)*v + dp)/(1+d) + 1/g, lh, ch, Gh)
    lo = np.zeros_like(d); hi = np.full_like(d, LT)
    for _ in range(64):
        mid = (lo + hi)/2; up = fa(mid) > fb(mid); lo = np.where(up, mid, lo); hi = np.where(up, hi, mid)
    G = (1 + d)*np.maximum(np.minimum(fa(lo), fb(lo)), np.minimum(fa(hi), fb(hi)))
    G0 = G[0]
    q = float(((G[2:] - G0 + 1e-12)/(DQ[1:-1]*G0)).max())
    q = max(q, float((G[1] - G0 + 1e-12)/(DQ[1]*G0)))
    return float(G0) + 1e-7, min(q, 1.0)
cache = {}
def T(lg):
    if lg not in cache: cache[lg] = Z.r3terms(lg, base)
    return cache[lg]
cand = []
for x0 in range(1, 17):
    K = 2**x0; cand.append((float(x0), (((K+1)**2)//4 - 1)/K**2, 1.0, 0))
t0 = time.time()
for lg in np.arange(6, 40.01, 0.25):
    for dd in (0.0, 0.5, 1.0, 2.0, 4.0):
        e, q = rowsolver(lg, lg + dd, T); cand.append((2*lg + dd, e, q, 1))
        if not TWO and not math.isinf(c):   # transposed grid (corner on the other child)
            e, q = rowsolver(lg + dd, lg, T); cand.append((2*lg + dd, e, q, 1))
print(f"c={cs} TWO={TWO} bases done {time.time()-t0:.0f}s", flush=True)
nb = int(math.log(XMAX/0.5)/math.log(ratio)) + 2
bx = np.full(nb, np.nan); be_ = np.full(nb, np.nan); bmu = np.zeros(nb); bq = np.ones(nb); bk = np.zeros(nb, int)
def put(xs, es, qs, ks):
    xs = np.ravel(xs); es = np.ravel(es); qs = np.ravel(qs); ks = np.ravel(ks)
    ok = (xs < XMAX) & (es > 0) & np.isfinite(es)
    xs, es, qs, ks = xs[ok], es[ok], qs[ok], ks[ok]
    if xs.size == 0: return
    mu = 1/(es*es*xs); b = np.clip((np.log(xs/0.5)/math.log(ratio)).astype(int), 0, nb-1)
    order = np.argsort(-mu, kind='stable'); b, mu, xs, es, qs, ks = b[order], mu[order], xs[order], es[order], qs[order], ks[order]
    _, first = np.unique(b, return_index=True)
    b, mu, xs, es, qs, ks = b[first], mu[first], xs[first], es[first], qs[first], ks[first]
    m = mu > bmu[b]; bmu[b[m]] = mu[m]; bx[b[m]] = xs[m]; be_[b[m]] = es[m]; bq[b[m]] = qs[m]; bk[b[m]] = ks[m]
put([c_[0] for c_ in cand], [c_[1] for c_ in cand], [c_[2] for c_ in cand], [c_[3] for c_ in cand])
lts = np.arange(1.0, 48.01, 1.0)
sdc = {}
def sdrow(x):
    if x not in sdc: sdc[x] = np.array([Z.sdS2(lt, x) for lt in lts])
    return sdc[x]
def zrow(i):
    x = bx[i]; K = 2.0**min(x, 60)
    return bq[i]*fcov(x, bk[i])*sdrow(x) + (1/(K-1) if x < 60 else 0.0)
def level(ea, eb, al, bet, ent):
    v, w = levH(ea, eb, al, bet)
    Ew = hsup(w*ea, (1-w)*eb); Sw = w*ea + (1-w)*eb; Ow = w*ea*al + (1-w)*eb*bet
    e2 = np.maximum(v, Ew + Ow) + ent
    q2 = np.minimum((Sw - Ew + Ow + ent)/e2, 1.0)
    return e2, q2
t0 = time.time()
for i in range(nb):
    if np.isnan(bx[i]): continue
    xa, ea = bx[i], be_[i]; za = zrow(i)
    js = [j for j in range(i + 1) if not np.isnan(bx[j]) and xa + bx[j] + 2 < XMAX]
    if not js: continue
    xb = bx[js][:, None, None]; eb = be_[js][:, None, None]
    zb = np.array([zrow(j) for j in js])[:, None, :]
    lg = xa + lts[None, :, None]; lh = xb + lts[None, None, :]
    al = 2.0**-lh + za[None, :, None]; bet = 2.0**-lg + zb
    ent = 2.0**-np.minimum(xa, xb)
    e2, q2 = level(ea, eb, al, bet, ent)
    if not TWO and not math.isinf(c):
        e3, q3 = level(eb, ea, bet, al, ent)
        m = e3 < e2; e2 = np.where(m, e3, e2); q2 = np.where(m, q3, q2)
    put(np.broadcast_to(lg + lh, e2.shape), e2, q2, np.ones(e2.shape, int))
ok = ~np.isnan(bx); X = bx[ok]; E = be_[ok]; Q = bq[ok]; KD = bk[ok]; o = np.argsort(X); X, E, Q, KD = X[o], E[o], Q[o], KD[o]
Em = np.minimum.accumulate(E)
tag = f"c{cs}{'_two' if TWO else ''}"
json.dump({'x': X.tolist(), 'e': E.tolist(), 'q': Q.tolist(), 'kind': KD.tolist()}, open(f'sol_{tag}.json', 'w'))
print(f"c={cs} time {time.time()-t0:.0f}s")
for Y in (16, 32, 64, 128, 256, 512, 1000, 2000, 4000, 8000, 16000):
    k = np.searchsorted(X, Y, side='right') - 1
    if k >= 0: print(f"  x<={Y}: e={Em[k]:.5f} 2e*sqrt(2x)={2*Em[k]*math.sqrt(2*Y):.3f} e*x/log2x={Em[k]*Y/math.log2(Y):.3f}")
for a in (2,):
    s = []
    for Y in (1000, 2000, 4000, 8000):
        ys = np.exp(np.linspace(math.log(Y), math.log(2*Y), 600)); k = np.searchsorted(X, ys, side='right') - 1
        s.append(f"[{Y},{2*Y}] {(2*Em[k]*np.sqrt(a*ys)).max():.3f}")
    print(f"  a={a}: " + "; ".join(s))
