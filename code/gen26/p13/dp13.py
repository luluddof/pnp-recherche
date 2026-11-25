# G26-P13: DP of the Lemma-U level model (P10's dp10.py) extended to EDGE-WEIGHTED children (structural change S2).
# Each solver carries (x, e, q, kind): cost on an edge-weighted clique <= e*(T + q*Dev), T = calibrated total,
#   Dev = Dbar (placement-averaged rectangle-diagonal differences) for kind 'rect' (levels, row-sweep bases),
#   Dev = D = 2 sum_e |lam_e - c| for kind 'plain' (clique sweep; q = 1 exactly).
# Parent grouping excess: z = q_child * f(kind_child, K) * sdS2(t, K) + 1/(K-1),
#   f = 1 (D-route, |l1-l2| <= |l1-c|+|l2-c|), or f = sqrt((1+eps_K)/2) (covariance lemma, rect children only).
# modes: pad (q=1, f=1: P10's model) | padiii (q=1, f) | split (q from splitting, f=1) | splitiii (q, f)
import math, sys, time, json
import numpy as np
import zlib10 as Z
base = sys.argv[1]; mode = sys.argv[2]; ratio = float(sys.argv[3]) if len(sys.argv) > 3 else 1.02
XMAX = 2.0**15
SPLIT = mode.startswith('split'); III = mode.endswith('iii')
def fcov(x, kind):
    if kind == 0 or not III: return 1.0
    K = 2.0**min(x, 60)
    if K < 8: return 1.0
    eps = (4*K - 2)/((K - 2)*(K - 3))
    return min(1.0, math.sqrt((1 + eps)/2))
def lev2(ea, eb, al, be):
    """max_th min(ea(cos th + al), eb(sin th + be)) (= P10's lev) and the minimax weight w."""
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
def sc(xr, l, cA, Gam):
    gg = 2.0**l; v = np.maximum(xr, 2*math.sqrt(cA))
    return xr*gg/(4*(gg-1)) + (v - xr)/4 + cA/v + Gam + 2/gg
THS = np.linspace(0, math.pi/2, 20001); THQ = np.linspace(0, math.pi/2, 4001)
DQ = np.concatenate([[0.0], np.geomspace(1e-5, 1e3, 160)])
def rowsolver(lg, lh, T):
    cg, Gg = T(lg); ch, Gh = T(lh); g, h = 2.0**lg, 2.0**lh
    e = float(np.minimum(sc(1/h + np.cos(THS), lg, cg, Gg), sc(1/g + np.sin(THS), lh, ch, Gh)).max()) + 1e-7
    if not SPLIT: return e, 1.0
    # G(d) = (1+d) max_th min(sc(xr_a), sc(xr_b)),  xr_a = ((1-d')cos th + d')/(1+d) + 1/h,  d' = min(d,1).
    # sc nondecreasing in xr; xr_a decreasing, xr_b increasing in th -> max-min at the crossing: bisection (60 steps).
    d = DQ; dp = np.minimum(d, 1.0)
    fa = lambda th: sc(((1-dp)*np.cos(th) + dp)/(1+d) + 1/h, lg, cg, Gg)
    fb = lambda th: sc(((1-dp)*np.sin(th) + dp)/(1+d) + 1/g, lh, ch, Gh)
    lo = np.zeros_like(d); hi = np.full_like(d, math.pi/2)
    for _ in range(60):
        mid = (lo + hi)/2; up = fa(mid) > fb(mid); lo = np.where(up, mid, lo); hi = np.where(up, hi, mid)
    G = (1 + d)*np.maximum(np.minimum(fa(lo), fb(lo)), np.minimum(fa(hi), fb(hi)))
    G0 = G[0]
    # G nondecreasing in d: chord slope on [d_i, d_{i+1}] <= (G(d_{i+1}) - G0)/d_i ; margin 1e-12 (bisection accuracy)
    q = float(((G[2:] - G0 + 1e-12)/(DQ[1:-1]*G0)).max())
    q = max(q, float((G[1] - G0 + 1e-12)/(DQ[1]*G0)))
    q = q*e/G0 if False else q
    return e, min(q, 1.0)     # per solver: split if better, else pad (q = 1 exactly)
cache = {}
def T(lg):
    if lg not in cache: cache[lg] = Z.r3terms(lg, base)
    return cache[lg]
cand = []   # (x, e, q, kind)
def ddc(K): return 0 if K <= 1 else (K//2)*(K - K//2) + ddc(K - K//2)
for x0 in range(1, 17):
    K = 2**x0; cand.append((float(x0), (ddc(K) if base == 'poly' else ((K+1)**2)//4 - 1)/K**2, 1.0, 0))
t0 = time.time()
for lg in (np.arange(6, 40.01, 0.25) if base != 'poly' else []):
    for dd in (0.0, 0.5, 1.0, 2.0, 4.0):
        e, q = rowsolver(lg, lg + dd, T); cand.append((2*lg + dd, e, q, 1))
print(f"bases done {time.time()-t0:.0f}s; q of row-sweep bases: min {min(c[2] for c in cand if c[3]==1):.4f} max {max(c[2] for c in cand if c[3]==1):.4f}", flush=True)
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
put([c[0] for c in cand], [c[1] for c in cand], [c[2] for c in cand], [c[3] for c in cand])
lts = np.arange(1.0, 48.01, 1.0)
sdc = {}
def sdrow(x):
    if x not in sdc: sdc[x] = np.array([Z.sdS2(lt, x) for lt in lts])
    return sdc[x]
def zrow(i):
    x = bx[i]; K = 2.0**min(x, 60)
    return bq[i]*fcov(x, bk[i])*sdrow(x) + (1/(K-1) if x < 60 else 0.0)
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
    v, w = lev2(ea, eb, al, bet)
    ent = 2.0**-np.minimum(xa, xb)
    Ew = np.hypot(w*ea, (1-w)*eb); Sw = w*ea + (1-w)*eb; Ow = w*ea*al + (1-w)*eb*bet
    e2 = np.maximum(v, Ew + Ow) + ent
    q2 = np.minimum((Sw - Ew + Ow + ent)/e2, 1.0) if SPLIT else np.ones_like(e2)
    put(np.broadcast_to(lg + lh, e2.shape), e2, q2, np.ones(e2.shape, int))
ok = ~np.isnan(bx); X = bx[ok]; E = be_[ok]; Q = bq[ok]; KD = bk[ok]; o = np.argsort(X); X, E, Q, KD = X[o], E[o], Q[o], KD[o]
Em = np.minimum.accumulate(E)
json.dump({'x': X.tolist(), 'e': E.tolist(), 'q': Q.tolist(), 'kind': KD.tolist()}, open(f'sol13_{base}_{mode}_{ratio}.json', 'w'))
print(f"base={base} mode={mode} ratio={ratio} time {time.time()-t0:.0f}s")
sel = (X > 100) & (KD == 1)
if sel.any(): print(f"  q of level solvers with x>100: min {Q[sel].min():.4f} median {np.median(Q[sel]):.4f} max {Q[sel].max():.4f}")
for a in (2, 4):
    s = []
    for Y in (1000, 2000, 4000, 8000):
        ys = np.exp(np.linspace(math.log(Y), math.log(2*Y), 600)); k = np.searchsorted(X, ys, side='right') - 1
        s.append(f"[{Y},{2*Y}] {(2*Em[k]*np.sqrt(a*ys)).max():.3f}")
    print(f"  a={a}: " + "; ".join(s))
