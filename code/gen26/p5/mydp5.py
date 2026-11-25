# Own DP over unbalanced solver trees (Lemma U) with R2 grouping and R3 row-sweep base. Asymptotic model (lam -> inf).
import math, sys, time, json
import numpy as np
import mylib as M
mode = sys.argv[1]; zname = sys.argv[2]; ratio = float(sys.argv[3]) if len(sys.argv) > 3 else 1.02
ltstep = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
XMAX = 2.0**15
zeta = {'hoeff': M.zeta_hoeff, 'R2': M.zeta_R2, 'F': M.zeta_F, 'Fnomcd': (lambda lt, x: M.zeta_F(lt, x, mcd=False)), 'floor': M.zeta_floor, 'lin2': M.zeta_lin2, 'randfloor': M.zeta_randfloor, 'S2': M.zeta_S2, 'zero': M.zeta_zero}[zname]
_zc = {}
def zrow(x):
    if x not in _zc: _zc[x] = np.asarray(zeta(lts, x), float)
    return _zc[x]
def p9_terms(lg):   # p9 Lemma RS as a (cA, Gam) pair: cA = 0, Gam = zs0 - 1/(4(g-1)) - 2/g + ... (keep zs total)
    g = 2.0**lg; return 0.0, math.sqrt(math.log(2*g*g+2*g+1)/(4*g))*(1+1/g)
def bases():
    S = []
    if mode.startswith('poly'):
        for x0 in range(1, 17): K = 2**x0; S.append((float(x0), M.d_dc(K)/K**2))
    else:
        for x0 in range(1, 17): K = 2**x0; S.append((float(x0), M.s_sweep(K)/K**2))
    if mode.startswith('row'):
        if mode.startswith('rowp9'): terms = p9_terms
        elif 'b7c' in mode: terms = lambda lg: M.r3_terms(lg, b7const=True)
        elif 'zero' in mode: terms = lambda lg: (0.0, 0.0)
        elif 'ser' in mode: terms = M.r3_terms_v2
        elif 'half' in mode: terms = lambda lg: (lambda r: (r[0], r[1]/2))(M.r3_terms(lg))   # R3 with Hoeffding-Azuma (interval) for Gamma
        else: terms = M.r3_terms
        cache = {}
        T = lambda lg: cache.setdefault(lg, terms(lg))
        for lg in np.arange(6, 40.01, 0.25):
            for d in (0.0, 0.5, 1.0, 2.0, 4.0):
                S.append((lg + lg + d, M.rowsweep_solver(lg, lg + d, T)))
    return S
nb = int(math.log(XMAX/0.5)/math.log(ratio)) + 2
bx = np.full(nb, np.inf); be = np.full(nb, np.inf); bmu = np.zeros(nb)
def binof(x): return np.clip((np.log(np.asarray(x)/0.5)/math.log(ratio)).astype(int), 0, nb-1)
def insert(xs, es):
    xs = np.ravel(xs); es = np.ravel(es); ok = (xs < XMAX) & np.isfinite(es) & (es > 0)
    xs, es = xs[ok], es[ok]
    if xs.size == 0: return
    mu = 1/(es*es*xs); b = binof(xs)
    o = np.lexsort((mu, b)); b, mu, xs, es = b[o], mu[o], xs[o], es[o]
    last = np.r_[b[1:] != b[:-1], True]; b, mu, xs, es = b[last], mu[last], xs[last], es[last]
    u = mu > bmu[b]; bmu[b[u]] = mu[u]; bx[b[u]] = xs[u]; be[b[u]] = es[u]
lts = np.arange(1.0, 48.01, ltstep)
t0 = time.time()
S = bases(); insert([s[0] for s in S], [s[1] for s in S])
LA, LB = np.meshgrid(lts, lts, indexing='ij')
for i in range(nb):
    if not np.isfinite(bx[i]): continue
    xa, ea = bx[i], be[i]; za = zrow(xa)[:, None]*np.ones_like(LA)
    XX = []; EE = []
    for j in range(i+1):
        if not np.isfinite(bx[j]): continue
        xb, eb = bx[j], be[j]
        if xa + xb + 2 >= XMAX: continue
        zb = zrow(xb)[None, :]*np.ones_like(LB)
        lg = xa + LA; lh = xb + LB
        e2 = M.level_U(ea, eb, np.power(2.0, -lh) + za, np.power(2.0, -lg) + zb) + 2.0**-min(xa, xb)
        XX.append((lg + lh).ravel()); EE.append(e2.ravel())
    if XX: insert(np.concatenate(XX), np.concatenate(EE))
fin = np.isfinite(bx); X = bx[fin]; E = be[fin]; o = np.argsort(X); X = X[o]; E = E[o]
Em = np.minimum.accumulate(E)
json.dump({'x': X.tolist(), 'e': E.tolist()}, open(f'solvers_{mode}_{zname}_{ratio}_{ltstep}.json', 'w'))
print(f"mode={mode} zeta={zname} ratio={ratio} ltstep={ltstep} time {time.time()-t0:.0f}s")
for a in (2, 4):
    out = []
    for Y in (1000, 2000, 4000, 8000):
        ys = np.exp(np.linspace(math.log(Y), math.log(2*Y), 600))
        k = np.searchsorted(X, ys, side='right') - 1
        v = 2*Em[k]*np.sqrt(a*ys); out.append(f"[{Y},{2*Y}] max {v.max():.3f} min {v.min():.3f}")
    print(f"  a={a}: " + "; ".join(out))
for xm in (100, 500, 2000, 8000, 16000):
    k = np.searchsorted(X, xm, side='right') - 1; print(f"  x<={xm}: e={Em[k]:.5f} mu={1/(Em[k]**2*X[k]):.4f}")
