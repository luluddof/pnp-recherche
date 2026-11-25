# Own threshold computation (structure of p2 thmA_explicit 'lasvegas' top, re-implemented), on a solver table (x, e).
import math, sys, json
import numpy as np
import mylib as M
def lp(n): return (math.pi*math.sqrt(2*n/3) - math.log(4*n*math.sqrt(3)))/math.log(2)   # HR leading term (p2 model)
def rho_bar_small(L2):
    lam = 2.0**L2; lo, hi = 1, 4
    while lp(hi) < lam: hi *= 2
    while lo < hi:
        mid = (lo+hi)//2
        if lp(mid) >= lam: hi = mid
        else: lo = mid+1
    h = lo; S = 2*math.ceil(math.log2(2*h+1) + lp(2*h) + 0.5)
    return S/(4*math.ceil(lam))
LN2 = math.log(2)
def top_cost(x, e, logm, r, topgroup='hoeff', depth=12):
    """exponent/n with solver (x = log2 K_top, e) at the top; derandomised Lemma 8' (R1) + Lemma 8 (ii) grouping."""
    ln_ = logm - math.log2(r); wn = 3*r - 1; lW = logm + math.log2(3 - 1/r)
    best = None
    for lt4 in range(2, 160):
        lt = lt4/2; lg = x + lt
        if lg + math.log2(x) > logm - 1: break
        z = float(M.zeta_hoeff(lt, x) if topgroup == 'hoeff' else M.zeta_R2(lt, x))
        lmu = 1 + lW - 2*lg
        lmu = min(lmu, 1000.0)
        ldc = 0.5*(math.log2(2*(4/3)*(math.log(4) + lg*LN2)) - (logm - lg))
        if ldc > 0: continue
        dc = 2.0**ldc
        for d in (0.01, 0.02, 0.05, 0.1, 0.2, 0.4):
            # Bernstein, 9 classes, (g^2/2) pairs, x9 extra (as p2): <= 1/4 ; derandomised by R1 (Phi(0) <= 1/2)
            if math.log(40.5) + 2*lg*LN2 - d*d*2.0**lmu/(18*(1+d/3)) > math.log(0.25): continue
            dev = z*(1+2*d)
            L2 = 1 + lW + math.log2(1+dev) - 2*x                    # log2 lam_top
            # rho_bar correction (conservative: every level at lam_top), p2's asymptotic form for log2 lam > 40
            rb = 1 + ((10*L2 + 30)/4*2.0**-L2)/M.RHO if L2 > 40 else (rho_bar_small(max(L2, 3))/M.RHO)
            c = 2*e*wn*(1+dev)*rb**depth + 3*r*(1+dc)*2.0**-x
            # lower order: grouping search g log K, |R| log|I| <= K^2 (2 L2 + 2), ceil(L): e K^2, poly per step
            for v in (lg + math.log2(x), 2*x + math.log2(2*max(L2, 1) + 2), math.log2(e) + 2*x, math.log2(100*(4*logm + 64))):
                c += 2.0**(v - ln_)
            if best is None or c < best[0]: best = (c, lt, d, L2)
    return best
def threshold(sol, r, tau, topgroup):
    X, E = sol
    def f(lm):
        best = None
        for x, e in zip(X, E):
            if 2*x + 2 > lm: break
            if 2*e*(3*r-1) > tau: continue
            c = top_cost(x, e, lm, r, topgroup)
            if c and (best is None or c[0] < best[0]): best = (c[0], x, e, c[1], c[2], c[3])
        return best
    lo, hi = 50.0, float(__import__("os").environ.get("THR_HI", "3e4"))
    b = f(hi)
    if b is None or b[0] > tau: return None
    for _ in range(30):
        mid = (lo+hi)/2; b = f(mid)
        if b is not None and b[0] <= tau: hi = mid
        else: lo = mid
    return hi, f(hi)
if __name__ == "__main__":
    fn = sys.argv[1]; topgroup = sys.argv[2] if len(sys.argv) > 2 else 'hoeff'
    if fn.endswith('.json'):
        d = json.load(open(fn)); X = np.array(d['x']); E = np.array(d['e'])
    else:   # p2's E2 chain for calibration: x0 = 8, zeta_l = 0.0625*0.5^l, t_l least with 2 z^2 t >= ln(2 t K^2)
        X = []; E = []; x = 8.0; K = 2**8; e = (((K+1)**2)//4 - 1)/K**2; X.append(x); E.append(e)
        for l in range(12):
            z = 0.0625*0.5**l; t = 1.0
            for _ in range(200): t = max(2.0, (math.log(2*t) + 2*x*LN2)/(2*z*z))
            t = math.ceil(t); lt = math.log2(t)
            e = e*(M.RHO + 2.0**-(x+lt) + float(M.zeta_hoeff(lt, x))) + 2.0**-x; x = 2*(x + lt); X.append(x); E.append(e)
        X = np.array(X); E = np.array(E)
    # keep Pareto-useful solvers only (decreasing e with x) and thin them
    o = np.argsort(X); X, E = X[o], E[o]; keep = []; em = np.inf
    for i in range(len(X)):
        if E[i] < em: keep.append(i); em = E[i]
    X, E = X[keep], E[keep]
    for r in (4.26, 1.0):
        row = []
        for tau, nm in ((1.0, '2^n'), (math.log2(1.6423), '1.6423^n'), (math.log2(1.30704), '1.30704^n')):
            res = threshold((X, E), r, tau, topgroup)
            if res is None: row.append(f"{nm}: none"); continue
            th, b = res
            row.append(f"{nm}: {th:.4g} (x_top={b[1]:.1f}, e={b[2]:.5f}, lt={b[3]}, d={b[4]}, log2lam={b[5]:.1f})")
        print(f"[{fn} top={topgroup}] r={r}: " + " | ".join(row), flush=True)
