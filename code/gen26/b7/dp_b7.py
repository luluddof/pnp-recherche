# G25-P9 levers (b) unbalanced grid levels and (d) sweep directly on the row network ("row sweep", level 1
# without grouping).  Asymptotic model (lam -> infinity: App. B frontier = quarter circle, rho_bar -> 1/sqrt2;
# per-level additive overheads kappa, c_ev ignored, as in p2/asym_c1.py).
#
# Solver = (x = log2 K, e).  Unbalanced grid level from solvers A (rows) and B (columns), grouping factors
# t_a, t_b (Lemma 8 (ii) with split zeta, z = zeta1+zeta2 = zsum(t, x)):
#   log g = x_a + log t_a,  log h = x_b + log t_b,  x' = log g + log h,
#   e' <= sqrt((w e_a)^2 + ((1-w) e_b)^2) + w e_a (1/h + z_a) + (1-w) e_b (1/g + z_b) + 1/min(K_a, K_b),
#   any w in [0,1]; we use w = e_b^2/(e_a^2+e_b^2).   (Lemma U in notes.md)
# Row-sweep solver (no grouping), g x h grid, sweep the row network on g nodes or the column network on h nodes:
#   e <= (1/4) sqrt(w^2+(1-w)^2) + w (1/(4h) + zs(g)) + (1-w) (1/(4g) + zs(h)),
#   zs(g) = sqrt(ln(2g^2+2g+1)/(4g)) (1+1/g) + 1/(4(g-1)) + 2/g                  (Lemma RS in notes.md)
# Clique-sweep base: e0 = (floor((K0+1)^2/4)-1)/K0^2 (p2 Lemma S).
import math, sys
import numpy as np
import b7lib as B
MODE = sys.argv[1] if len(sys.argv) > 1 else 'new'
LN2 = math.log(2)

def zsum_np(lt, x):
    if MODE in ('new','poly'): return B.zeta_new(lt, x)   # split zeta of p2 (thmA_explicit.zsum with split=True); lt = log2 t, x = log2 K
    t = 2.0**lt; lnK = x*LN2
    return np.sqrt((math.log(2.0001) + lnK)/(2*t)) + np.sqrt((np.log(2.0001*t) + 2*lnK)/(2*t))

def zs_sweep(lg):
    if MODE == 'new': return B.zs_R3(float(lg))
    g = 2.0**lg
    return np.sqrt(np.log(2*g*g + 2*g + 1)/(4*g))*(1 + 1/g) + 1/(4*(g-1)) + 2/g

def base_solvers(rowsweep=True, xs=range(2, 17)):
    S = []
    for x0 in xs:
        K = 2**x0
        if MODE == 'poly':
            import chain_b7 as C; S.append((float(x0), C.dC(K)/K**2, 'clique-DC'))
        else: S.append((float(x0), (((K+1)**2)//4 - 1)/K**2, 'clique-sweep'))
    if rowsweep:
        for lg in np.arange(8, 40.01, 0.25):
            for d in (0.0, 1.0, 2.0, 4.0):
                lh = lg + d
                g, h = 2.0**lg, 2.0**lh
                a1 = 1/(4*h) + float(zs_sweep(lg)); b1 = 1/(4*g) + float(zs_sweep(lh))
                best = min(0.25*math.sqrt(w*w + (1-w)**2) + w*a1 + (1-w)*b1 for w in np.linspace(0, 1, 201))
                S.append((lg + lh, best, 'row-sweep'))
    return S

def run(rowsweep, unbalanced, a_list=(2, 4), xmax=2**15, ratio=1.02, lts=np.arange(4, 46, 1.0), verbose=True):
    nb = int(math.log(xmax/2)/math.log(ratio)) + 2
    bx = np.full(nb, np.inf); be = np.full(nb, np.inf); bmu = np.zeros(nb)
    def binof(x): return np.clip((np.log(x/2)/math.log(ratio)).astype(int), 0, nb-1)
    def insert(xs, es):
        xs = np.asarray(xs, float).ravel(); es = np.asarray(es, float).ravel()
        ok = xs < xmax
        xs, es = xs[ok], es[ok]
        if xs.size == 0: return
        mu = 1/(es*es*xs); b = binof(xs)
        o = np.lexsort((mu, b)); b, mu, xs, es = b[o], mu[o], xs[o], es[o]
        last = np.r_[b[1:] != b[:-1], True]
        b, mu, xs, es = b[last], mu[last], xs[last], es[last]
        upd = mu > bmu[b]
        bmu[b[upd]] = mu[upd]; bx[b[upd]] = xs[upd]; be[b[upd]] = es[upd]
    S = base_solvers(rowsweep)
    insert([s[0] for s in S], [s[1] for s in S])
    LTa, LTb = np.meshgrid(lts, lts, indexing='ij')
    for i in range(nb):
        if not np.isfinite(bx[i]): continue
        xa, ea = bx[i], be[i]
        js = [j for j in range(i+1) if np.isfinite(bx[j])] if unbalanced else [i]
        XX = []; EE = []
        for j in js:
            xb, eb = bx[j], be[j]
            za = zsum_np(LTa, xa); zb = zsum_np(LTb, xb)
            lg = xa + LTa; lh = xb + LTb
            w = eb*eb/(ea*ea + eb*eb)
            e2 = (np.sqrt((w*ea)**2 + ((1-w)*eb)**2) + w*ea*(2.0**-lh + za) + (1-w)*eb*(2.0**-lg + zb)
                  + 2.0**-min(xa, xb))
            XX.append((lg + lh).ravel()); EE.append(e2.ravel())
        insert(np.concatenate(XX), np.concatenate(EE))
    # c1 over ranges: best e among solvers with a x <= log m ; c1(m) = 2 e sqrt(log m)
    fin = np.isfinite(bx); X = bx[fin]; E = be[fin]
    o = np.argsort(X); X = X[o]; E = np.minimum.accumulate(E[o])
    res = {}
    for a in a_list:
        for Y in (2000, 4000, 8000):
            ys = np.exp(np.linspace(math.log(Y), math.log(2*Y), 400))
            vals = []
            for y in ys:
                k = np.searchsorted(X, y, side='right') - 1
                vals.append(2*E[k]*math.sqrt(a*y))
            res[(a, Y)] = (max(vals), min(vals))
    if verbose:
        mus = 1/(E*E*X)
        for xm in (500, 1000, 2000, 4000, 8000, 16000):
            k = np.searchsorted(X, xm, side='right') - 1
            print(f"   x<={xm:6d}: best e={E[k]:.5f} at x={X[k]:.1f}, mu={1/(E[k]**2*X[k]):.4f}")
    return res

if __name__ == "__main__":
    runs = [("POLY SPACE unbalanced, clique D&C base", False, True)] if MODE == 'poly' else [("unbalanced levels + row-sweep (MODE=%s)" % MODE, True, True),
                         ("unbalanced levels, clique-sweep base", False, True)]
    for name, rs, ub in runs:
        print(name, flush=True)
        r = run(rs, ub)
        for k, v in sorted(r.items()):
            print(f"   a={k[0]} log m/a in [{k[1]},{2*k[1]}]:  max c1 = {v[0]:.3f}   min = {v[1]:.3f}")
        sys.stdout.flush()
