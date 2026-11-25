# G25-P8: independent re-derivation of the explicit Theorem A' recursion, c1 and thresholds.
# Own code; only formulas re-derived from the paper's App. A text (Lemma 8, Thm A (40)-(45)).
import math, sys, itertools, random
LN2 = math.log(2); RHO = 1 / math.sqrt(2)

def zeta12(t, x):
    """split Hoeffding (Lemma 8 (ii), t = g/K, K = 2^x): K e^{-2t z1^2} + t K^2 e^{-2t z2^2} < 1."""
    z1 = math.sqrt((math.log(2.0001) + x * LN2) / (2 * t))
    z2 = math.sqrt((math.log(2.0001 * t) + 2 * x * LN2) / (2 * t))
    return z1, z2

def zeta_paper(t, x):
    """paper's single zeta: 2 z^2 g >= K ln(2gK), i.e. z = sqrt(ln(2 t K^2)/(2t)); deviation 2z."""
    return math.sqrt((math.log(2 * t) + 2 * x * LN2) / (2 * t))

def step(e, x, lt, rb=RHO, split=True, entry=1, extra=0.0):
    t = 2.0 ** lt
    if split: z1, z2 = zeta12(t, x); dev = z1 + z2
    else: dev = 2 * zeta_paper(t, x)
    ig = 2.0 ** (-(x + lt)) if x + lt < 1000 else 0.0
    iK = 2.0 ** (-x) if x < 1000 else 0.0
    e2 = e * (rb + ig + dev) + entry * iK + extra * ig ** 4
    return e2, 2 * (x + lt)

def run(x0, lts, base='sweep', **kw):
    K0 = 2 ** x0
    e = ((K0 + 1) ** 2 // 4 - 1) / K0 ** 2 if base == 'sweep' else (K0 - 1) / (2 * K0)
    xs, es = [x0], [e]; x = x0
    for lt in lts:
        e, x = step(e, x, lt, **kw); xs.append(x); es.append(e)
    return xs, es

# ---------- paper as written (40) ----------
def paper40(x0=18, L=30):
    xs, es = [x0], [0.5]; x, e = x0, 0.5
    for l in range(L):
        z = 2.0 ** (-l - 5)
        e = (RHO + 3 * z) * e + (3 * 2.0 ** (-x) if x < 1000 else 0)
        x = 2 * (x + 3 * l + x0)          # g_l = 2^{3l+x0} K_l, K_{l+1} = g_l^2
        xs.append(x); es.append(e)
    return xs, es

# ---------- top levels: exponent / n at given log2 m (r = m/n, W = 3m - n) ----------
def top_E1(x, e, logm, r):
    """Lemma T4 (deterministic Chebyshev, refined variance). returns best exponent/n or None"""
    n = 2.0 ** (logm - math.log2(r)); m = 2.0 ** logm; W = 3 * m - n; K = 2.0 ** x
    best = None
    for k in range(0, 400):
        lg = 4 * x + k / 4
        if lg > logm - 2: break
        g = 2.0 ** lg; P = 6 * math.ceil(m / g) if m / g < 1e300 else 6 * m / g
        # pairs: P (1+2/K) K^4 / (2 zp^2 * 2W) <= 1/4
        zp2 = 4 * P * (1 + 2 / K) * K ** 4 / (4 * W) if x < 200 else 2.0 ** (math.log2(P) + 4 * x - math.log2(W))
        zp = math.sqrt(zp2)
        if zp >= 1: continue
        dp = math.sqrt(3 * K * math.log(4 * K) / g) if x < 200 else 2.0 ** (0.5 * (math.log2(3 * (math.log(4) + x * LN2)) + x - lg))
        if dp > 1: continue
        c = (e * (2 * W * (1 + zp)) + 3 * (1 + dp) * (m + g) / K + g * x) / n
        if best is None or c < best: best = c
    return best

def top_logs(x, e, logm, r, variant):
    # numerically stable version for huge numbers: everything in log2, divided by n
    ln_ = logm - math.log2(r)
    wn = 3 * r - 1                       # W / n
    if variant == 'paper9':
        if logm < 9 * x: return None
        # Lemma 8 (i): L = (U + gP/K)/K^2, entry 2^{2gP/K}, gP <= 6(m+g), g = 3K^8, U <= 2W
        gm = 3 * 2.0 ** (8 * x - logm)
        return 2 * e * wn + e * 6 * r * (1 + gm) * 2.0 ** -x + 12 * r * (1 + gm) * 2.0 ** -x + 3 * 2.0 ** (8 * x + math.log2(x) - ln_)
    if variant == 'E1':
        best = None
        for k in range(0, 600):
            lg = 4 * x + k / 4
            if lg + math.log2(x) > ln_: break
            # zp^2 = P (1+2/K) K^4 / W, P = 6 m/g (1+g/m)
            lzp2 = math.log2(6) + logm - lg + math.log2(1 + 2.0 ** (lg - logm)) + math.log2(1 + 2 * 2.0 ** -x) + 4 * x - (ln_ + math.log2(wn))
            if lzp2 >= 0: continue
            zp = 2.0 ** (lzp2 / 2)
            ldp2 = math.log2(3 * (math.log(4) + x * LN2)) + x - lg
            if ldp2 > 0: continue
            dp = 2.0 ** (ldp2 / 2)
            c = 2 * e * wn * (1 + zp) + 3 * r * (1 + dp) * (1 + 2.0 ** (lg - logm)) * 2.0 ** -x + 2.0 ** (lg + math.log2(x) - ln_)
            if best is None or c < best: best = c
        return best
    if variant == 'E2':
        best = None
        lW = ln_ + math.log2(wn)
        for k in range(2, 600):
            lt = k / 4; lg = x + lt
            if lg + math.log2(x) > ln_: break
            t = 2.0 ** lt
            z1, z2 = zeta12(t, x)
            lmu = 1 + lW - 2 * lg                # mu_bar = 2W/g^2
            if lmu > 900: mu = float('inf')
            else: mu = 2.0 ** lmu
            for d in [0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8]:
                # (9/2) g^2 exp(-d^2 mu/(18(1+d/3))) <= 1/4
                if math.log(4.5) + 2 * lg * LN2 - d * d * mu / (18 * (1 + d / 3)) > math.log(0.25): continue
                # g exp(-dc^2 (m/g)/(2(1+dc/3))) <= 1/4, dc <= 1:  dc^2 >= 2(4/3) ln(4g)/(m/g)
                ldc2 = math.log2(8 / 3 * (math.log(4) + lg * LN2)) - (logm - lg)
                if ldc2 > 0: continue
                dc = 2.0 ** (ldc2 / 2)
                gP_over_n = max((1 + 2 * d) * 2 * wn, 3 * (1 + dc) * r)
                main = e * (2 * wn + (z1 + z2) * gP_over_n)
                c = main + 3 * r * (1 + dc) * 2.0 ** -x + 2.0 ** (lg + math.log2(x) - ln_) + 2.0 ** -ln_
                if best is None or c < best: best = c
        return best

def threshold(x, e, r, tau, variant, lo=10.0, hi=3e6):
    f = lambda lm: top_logs(x, e, lm, r, variant)
    c = f(hi)
    if c is None or c > tau: return None
    for _ in range(50):
        mid = (lo + hi) / 2; c = f(mid)
        if c is not None and c <= tau: hi = mid
        else: lo = mid
    return hi

def best_threshold(r, tau, variant, base='sweep', split=True, entry=1, x0s=range(4, 16), maxlev=10, iters=300, seed=1):
    """optimise schedules: x0 and per-level log2 t_l by (i) geometric zeta seeds, (ii) coordinate descent."""
    a = {'paper9': 9, 'E1': 4, 'E2': 2}[variant]
    estar = tau / (2 * (3 * r - 1))
    def score(x0, lts):
        xs, es = run(x0, lts, base=base, split=split, entry=entry)
        if es[-1] >= estar: return None, xs, es
        th = threshold(xs[-1], es[-1], r, tau, variant)
        return th, xs, es
    best = None
    rng = random.Random(seed)
    for x0 in x0s:
        for L in range(0, maxlev + 1):
            # seeds: constant step in log2 t
            for lt0 in [8, 10, 12, 14, 16, 18, 20, 22]:
                for inc in [0, 1, 2, 3]:
                    lts = [lt0 + inc * l for l in range(L)]
                    th, xs, es = score(x0, lts)
                    if th is not None and (best is None or th < best[0]): best = (th, x0, lts, xs, es)
    # coordinate descent from best
    th, x0, lts, xs, es = best
    improved = True
    while improved:
        improved = False
        for i in range(len(lts)):
            for dlt in [-2, -1, -0.5, -0.25, 0.25, 0.5, 1, 2]:
                l2 = list(lts); l2[i] = max(1, l2[i] + dlt)
                t2, xs2, es2 = score(x0, l2)
                if t2 is not None and t2 < th - 1e-9:
                    th, lts, xs, es = t2, l2, xs2, es2; improved = True
        for dx in [-1, 1]:
            t2, xs2, es2 = score(x0 + dx, lts)
            if t2 is not None and t2 < th - 1e-9:
                th, x0, xs, es = t2, x0 + dx, xs2, es2; improved = True
    return th, x0, lts, xs, es

def c1_limsup(base, split, entry, a, x0, ztarget, L=40):
    """schedule: zeta-ish target per level: choose lt so that z1+z2 ~ ztarget(l) ; returns 2 e_l sqrt(a x_{l+1})."""
    xs, es = [x0], [((2**x0 + 1) ** 2 // 4 - 1) / 4**x0 if base == 'sweep' else (2**x0 - 1) / 2**(x0 + 1)]
    x, e = xs[0], es[0]
    for l in range(L):
        zt = ztarget(l)
        # smallest lt (step 1/64) with z1+z2 <= zt
        lt = 1.0
        while True:
            z1, z2 = zeta12(2.0 ** lt, x)
            if z1 + z2 <= zt: break
            lt += 1 / 64
        e, x = step(e, x, lt, split=split, entry=entry)
        xs.append(x); es.append(e)
    return [2 * es[l] * math.sqrt(a * xs[l + 1]) for l in range(L)], xs, es

if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if what in ('paper', 'all'):
        xs, es = paper40()
        print("paper (40), x0=18: l, log2K, e:", [(l, xs[l], round(es[l], 5)) for l in range(12)])
        v = [2 * es[l] * math.sqrt(9 * xs[l + 1]) for l in range(29)]
        print("paper c1 seq 2 e_l sqrt(9 x_{l+1}) l=5,10,20,28:", [round(v[l], 3) for l in (5, 10, 20, 28)])
        vl = [2 * es[l] * math.sqrt(9 * xs[l]) for l in range(1, 29)]
        print("paper liminf seq 2 e_l sqrt(9 x_l) l=10,20,28:", [round(vl[l - 1], 3) for l in (10, 20, 28)], " ratio", round(vl[27] / v[27], 4))
        for r in (4.26, 1.0):
            # paper row A: main term only, W = 3m: need 6 r e_l <= 1, log m >= 9 x_l
            l = next(l for l in range(30) if 6 * r * es[l] <= 1)
            print(f"  row A r={r}: first l with 6 r e_l<=1: l={l}, e={es[l]:.5f}, 9 x_l = {9*xs[l]}")
            lB = next(l for l in range(30) if 2 * (3 * r - 1) * es[l] <= 1)
            print(f"  row B r={r}: W=3m-n: l={lB}, 9 x_l = {9*xs[lB]}")
