# Explicit, re-optimised Theorem A (Claude Oct-2026 preprint, Appendix A) -- thresholds with all terms counted.
# Exponent = log2(#operations) for #3-SAT of a 3-CNF with m = r n clauses, every variable occurring,
# W = 3m - n equality tests ("wires").  See notes.md for the lemmas behind each formula.
import math, sys, functools
LN2 = math.log(2); RHO = 1/math.sqrt(2)

# ---------------- Appendix B pieces: explicit row+column bit count S(lam) and |I| --------------------------
# c_s^row = ceil(log2 sum_{i=0}^{2h} p(s+i)), c_s^col = ceil(log2 sum_i p(2h-s-i)), h least with p(h) >= 2^lam.
# Paper bound used: R_s^row <= (2h+1) p(2h+s), R_s^col <= (2h+1) p(2h-s); log p concave => max at s = 0.
def lp_asym(n):   # log2 p(n), Hardy-Ramanujan leading term (relative error O(1/sqrt n))
    return (math.pi*math.sqrt(2*n/3) - math.log(4*n*math.sqrt(3)))/LN2
def h_of(lam):
    lo, hi = 1, 4
    while lp_asym(hi) < lam: hi *= 2
    while lo < hi:
        mid = (lo+hi)//2
        if lp_asym(mid) >= lam: hi = mid
        else: lo = mid+1
    return lo
def S_bound(lam):  # upper bound on max_s (c_s^row + c_s^col), +1 bit safety for the asymptotic formula
    h = h_of(lam)
    return 2*math.ceil(math.log2(2*h+1) + lp_asym(2*h) + 0.5)
def rho_bar_l(L2):  # L2 = log2 lam ; per-rectangle cost / (4 lam)
    if L2 > 40: return RHO + (10*L2 + 30)/4 * 2.0**-L2     # asymptotic regime (excess ~ c log lam)
    lam = 2.0**L2; return S_bound(math.ceil(lam))/(4*lam)
def rho_bar(lam): return rho_bar_l(math.log2(lam))
def logI_l(L2):
    if L2 > 40: return 2*L2 - 1.77 + 0.01
    return math.log2(4*h_of(math.ceil(2.0**L2))+1)
KAPPA = 100          # exponent of the polynomial cost per elementary step (additive per level) -- generous
CEVAL = 2*10         # pieces computed within Delta^10 = 2^{20 lam} operations (Def. 7(c)); additive per level

# ---------------- Lemma 8 (ii) balancing: exact Hoeffding, optionally split zeta ---------------------------
def zsum(t, x, split):   # x = log2 K, g = tK. Events: K group totals (range P), gK node-group sums (range P/g):
    # failure <= K e^{-2t z1^2} + tK^2 e^{-2t z2^2} < 1 (split) ; returns z1+z2
    lnK = x*LN2
    if split:
        z1 = math.sqrt((math.log(2.0001) + lnK)/(2*t)); z2 = math.sqrt((math.log(2.0001*t) + 2*lnK)/(2*t))
        return z1+z2
    z = math.sqrt((math.log(2*t) + 2*lnK)/(2*t))   # paper's condition 2 z^2 g >= K ln(2gK), z1 = z2 = z
    return 2*z

def t_for(z, x):  # least t with 2 z^2 t >= ln(2 t K^2) (paper's condition)
    t = 1.0
    for _ in range(200): t = max(2.0, (math.log(2*t) + 2*x*LN2)/(2*z*z))
    return math.ceil(t)

# ---------------- the recursion ---------------------------------------------------------------------------
def levels(x0, ts, base, split, entry1, rb):
    """x_l = log2 K_l, e_l.  ts[l] = g_l/K_l.  rb = rho_bar used at every level (valid: lam only grows downward).
       base: 'brute' e0 = (K0-1)/(2K0) ; 'sweep' e0 = (floor((K0+1)^2/4)-1)/K0^2.
       e_{l+1} = e_l (rb + 1/g + z1 + z2) + (1 or 2)/K + CEVAL/g^4 ."""
    K0 = 2**x0
    e = (K0-1)/(2*K0) if base == 'brute' else (((K0+1)**2)//4 - 1)/K0**2
    out = [(x0, e, None, None)]; x = x0
    for l, t in enumerate(ts):
        zs = zsum(t, x, split); g = t*2.0**x if x < 1000 else float('inf')
        inv_g = 1/g if x < 1000 else 0.0; invK = 2.0**-x if x < 1000 else 0.0
        e = e*(rb + inv_g + zs) + (1 if entry1 else 2)*invK + (CEVAL*inv_g**4)
        x = 2*(x + math.log2(t))
        out.append((x, e, t, zs))
    return out

def p2(v): return 2.0**v if v < 1000 else float('inf')

def overhead_per_n(lv, L2top, rb, ln_):
    """Sum over levels of the lower-order terms divided by n (log domain), going down from the top clique
       (K_L, lam_top), L2top = log2 lam_top, ln_ = log2 n.  Per level j>=1 (clique K_j, lam_j, from level j-1
       with g = t K_{j-1}): |R| log2|I(lam_j)| (choices iota) + g log2 K_{j-1} (grouping search)
       + e_{j-1} K_{j-1}^2 (ceil of L) + CEVAL lam_j (piece evaluation) + KAPPA log2 Xi_j ;
       level 0: log2 K0 + KAPPA log2 Xi_0 (sweep: K0 steps, poly cost per term)."""
    tot = 0.0; L2 = L2top
    for j in range(len(lv)-1, 0, -1):
        x, e, t, zs = lv[j]; xp, ep = lv[j-1][0], lv[j-1][1]
        lg = xp + math.log2(t)
        logXi = 4*(x + L2) + 64
        terms = [4*lg - 2 + math.log2(logI_l(L2)), lg + math.log2(xp), math.log2(ep) + 2*xp,
                 math.log2(CEVAL) + L2, math.log2(KAPPA*logXi)]
        tot += sum(p2(v - ln_) for v in terms)
        g_inv = 2.0**-lg
        L2 = math.log2(rb + g_inv + zs) + 4*lg + L2 - 2*xp
    tot += p2(math.log2(lv[0][0] + KAPPA*(4*(lv[0][0] + L2) + 64)) - ln_)
    return tot

# ---------------- top level variants ------------------------------------------------------------------------
def top_cost(variant, lv, logm, r):
    """exponent/n for the top level using the last entry of lv, or None if infeasible at this m."""
    x, e = lv[-1][0], lv[-1][1]
    ln_ = logm - math.log2(r); wn = 3*r - 1; lW = logm + math.log2(3 - 1/r)
    invK = p2(-x); base_ov = p2(math.log2(e) + 2*x - ln_) + p2(math.log2(KAPPA*(4*logm + 64)) - ln_)
    if variant == 'paper9':     # Lemma 8 (i) as in the paper: g = 3K^8, zeta = 1/(2K), K <= m^{1/9}
        if logm < 9*x: return None
        lg = math.log2(3) + 8*x; gm = p2(lg - logm)
        L2 = 1 + lW - 2*x
        c = 2*e*wn + e*6*r*(1+gm)*invK + 6*r*(1+gm)*invK + p2(lg + math.log2(x) - ln_) + base_ov
        return c + overhead_per_n(lv, L2, rho_bar_l(L2), ln_)
    if variant == 'cheb4':      # refined Chebyshev, deterministic search over all K^g groupings (Lemma 8 (i'))
        best = None
        for d4 in range(0, 120):
            lg = 4*x + d4/2
            if lg + math.log2(x) > logm - 1: break
            gm = p2(lg - logm)
            # P = 6 ceil(m/g) <= 6 (m/g)(1 + g/m);  Var/t^2 union over pairs <= P (1+2/K) K^4 / (4 zp^2 W) <= 1/4
            lzp = 0.5*(math.log2(6) + logm - lg + math.log2(1+gm) + math.log2(1+2*invK) + 4*x - lW)
            if lzp > -1: continue
            zp = 2.0**lzp
            ldp = 0.5*(math.log2(3*math.log(4)+3*x*LN2) + x - lg)   # Chernoff for group sizes, union <= 1/4
            if ldp > 0: continue
            dp = 2.0**ldp
            L2 = 1 + lW + math.log2(1+zp) - 2*x
            c = 2*e*wn*(1+zp) + 3*r*(1+dp)*(1+gm)*invK + p2(lg + math.log2(x) - ln_) + base_ov
            c += overhead_per_n(lv, L2, rho_bar_l(L2), ln_)
            if best is None or c < best: best = c
        return best
    if variant == 'lasvegas':   # random clustering of clauses + Lemma 8 (ii) at the top (Lemma 8')
        best = None
        for lt4 in range(2, 140):
            lt = lt4/2; t = 2.0**lt; lg = x + lt
            if lg + math.log2(x) > logm - 1: break
            zs = zsum(t, x, True)
            lmu = 1 + lW - 2*lg
            for d in [0.02, 0.05, 0.1, 0.2, 0.4]:
                # P(some pair load > (1+2d) 2W/g^2) <= (g^2/2) * 9 exp(-d^2 mu/(18(1+d/3))) <= 1/4
                lnpc = math.log(40.5) + 2*lg*LN2 - (d*d*p2(lmu)/(18*(1+d/3)) if lmu < 1000 else float('inf'))
                if lnpc > math.log(0.25): continue
                # cluster sizes <= (1+dc) m/g : g exp(-dc^2 (m/g)/(2(1+dc/3))) <= 1/4 (dc <= 1)
                ldc = 0.5*(math.log2(2*(4/3)*(math.log(4) + lg*LN2)) - (logm - lg))
                if ldc > 0: continue
                dc = 2.0**ldc
                dev = zs*(1+2*d)
                L2 = 1 + lW + math.log2(1+dev) - 2*x
                c = 2*e*wn*(1+dev) + 3*r*(1+dc)*invK + p2(lg + math.log2(x) - ln_) + base_ov
                c += overhead_per_n(lv, L2, rho_bar_l(L2), ln_)
                if best is None or c < best: best = c
        return best
    raise ValueError

def threshold_for(variant, lv, r, tau, lo=8, hi=2e6):
    """least log2 m (bisection, cost/n is decreasing in m) with exponent/n <= tau at the top level of lv."""
    f = lambda lm: top_cost(variant, lv, lm, r)
    c = f(hi)
    if c is None or c > tau: return None
    for _ in range(36):
        mid = (lo+hi)/2; c = f(mid)
        if c is not None and c <= tau: hi = mid
        else: lo = mid
    return hi

def sched_ts(x0, nlev, z0, q, ratio_mode):
    ts = []; x = x0
    for l in range(nlev):
        z = z0*q**l
        t = 2**(3*l+18) if ratio_mode == 'paper' else t_for(z, x)
        ts.append(t); x = 2*(x + math.log2(t))
    return ts

def optimise(variant, r, tau, base, split, entry1, ratio_mode='hoeff', x0s=range(3, 21), ntop=12):
    """Schedules: x0, number of levels, zeta_l = z0 q^l (t_l least with the paper's Hoeffding condition).
       Pre-filter by the main term (with the deviation slack of the top), keep the ntop schedules with the
       smallest log2 K_top, run the exact top-level cost on those."""
    rb0 = RHO + 1e-4
    cands = []
    for x0 in x0s:
        for nlev in range(0, 12):
            zgrid = [(2.0**-k, q) for k in range(2, 11) for q in (0.5, 0.6, 0.7, 0.8, 0.9)] if nlev > 0 else [(0, 0)]
            for z0, q in zgrid:
                ts = sched_ts(x0, nlev, z0, q, ratio_mode) if nlev else []
                lv = levels(x0, ts, base, split, entry1, rb0)
                if 2*lv[-1][1]*(3*r-1) > tau: continue
                cands.append((lv[-1][0] + 40*(2*lv[-1][1]*(3*r-1) > 0.97*tau), x0, nlev, z0, q, ts))
    cands.sort(key=lambda c: c[0])
    best = None
    for _, x0, nlev, z0, q, ts in cands[:ntop]:
        lv = levels(x0, ts, base, split, entry1, rb0)
        th = threshold_for(variant, lv, r, tau)
        if th is None: continue
        L2 = 1 + th + math.log2(3-1/r) - 2*lv[-1][0]
        lv2 = levels(x0, ts, base, split, entry1, rho_bar_l(max(L2, 6)))
        th2 = threshold_for(variant, lv2, r, tau)
        if th2 is None: continue
        if best is None or th2 < best[0]:
            best = (th2, dict(x0=x0, levels=nlev, z0=z0, q=q, log2K=[round(a[0], 1) for a in lv2],
                              e=[round(a[1], 5) for a in lv2], ts=ts, lam_top=L2))
    return best

if __name__ == "__main__":
    # monotonicity / size of rho_bar
    print("rho_bar(lam)/(1/sqrt2):", [(l, round(rho_bar(l)/RHO, 5)) for l in [64, 256, 1024, 4096, 2**14, 2**16, 2**20, 2**30]])
    rows = [
      ("C  I3: exact Hoeffding recursion, paper top (K<=m^1/9), brute base", 'paper9', 'brute', False, False),
      ("C' I3+: split zeta, entry 2^{N/K} (case ii)",                        'paper9', 'brute', True,  True),
      ("D  +I2: sweep base e0=1/4+1/(2K0)",                                   'paper9', 'sweep', True,  True),
      ("E1 +I4a: refined Chebyshev top (deterministic, a~4)",                 'cheb4',  'sweep', True,  True),
      ("E2 +I4: random clustering top, Las Vegas (a~2)",                      'lasvegas','sweep', True, True),
    ]
    which = sys.argv[1:] or ['4.26', '1']
    for rs in which:
        r = float(rs)
        print(f"\n=== density r = {r}  (W/n = {3*r-1:.2f})")
        for tau, tn in [(1.0, 'brute 2^n'), (math.log2(1.6423), '#3SAT 1.6423^n'), (math.log2(1.30704), 'PPSZ 1.30704^n')]:
            print(f" target {tn} (tau={tau:.4f})", flush=True)
            for name, var, base, split, e1 in rows:
                b = optimise(var, r, tau, base, split, e1)
                if b is None: print(f"   {name:62s} none"); continue
                d = b[1]
                print(f"   {name:62s} log2 m* = {b[0]:9.4g}  levels={d['levels']} x0={d['x0']} z0={d['z0']} q={d['q']} "
                      f"log2K={d['log2K']} e={d['e']} log2lam_top~{d['lam_top']:.3g}", flush=True)
