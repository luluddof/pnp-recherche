# Asymptotic cheap-ring corner curves for P_3 (normalised: rank = N^{3u}, ring = N^{3v}, N = 2^k):
#  (a) composite of two Majorana Bell windows (rank = N (p1+p2), ring = q1 q2), (b) native symmetric P(D),
#  (c) Theorem S'_3 lower bound (site-symmetric, d=2): rank >= N^3 Pr[L1 > max(L2,L3) + u(q)].
import math, sys
from math import lgamma, log
k = int(sys.argv[1]); L2 = k*log(2)
def lb(n, j): return lgamma(n+1) - lgamma(j+1) - lgamma(n-j+1)
def lse(xs):
    xs = [x for x in xs if x > -1e300]
    if not xs: return -1e301
    m = max(xs); return m + math.log(sum(math.exp(x-m) for x in xs))
def windows(lp, lc):
    import numpy as np
    lp = np.array(lp); K = len(lp) - 1
    pre = np.logaddexp.accumulate(lp)                 # log P[i <= s]
    suf = np.logaddexp.accumulate(lp[::-1])[::-1]      # log P[i >= s]
    clc = np.logaddexp.accumulate(np.array(lc))
    out = []
    for d in range(0, K):
        s = np.arange(0, K - d + 1)
        left = np.where(s >= 1, pre[np.maximum(s-1, 0)], -np.inf)
        right = np.where(s + d + 1 <= K, suf[np.minimum(s+d+1, K)], -np.inf)
        e = np.logaddexp(left, right)
        out.append((float(clc[min(d, len(clc)-1)]), float(np.min(e))))
    return out
# (a) Majorana: u in 0..k, P(u) = C(2k+1,2u+1)/4^k ; cost C(2k+1,j)
lpM = [lb(2*k+1, 2*u+1) - k*log(4) for u in range(k+1)]
lcM = [lb(2*k+1, j) for j in range(k+1)]
cM = windows(lpM, lcM)
compo = [(2*q/(3*L2), -(log(2)+e)/(3*L2)) for q, e in cM]          # (v, 1-u)
# (b) native: D = a - b, (a,b,c) ~ Mult(k; 1/4,1/4,1/2)
lpD = {}
for a in range(k+1):
    for b in range(k+1-a):
        c = k - a - b
        x = lgamma(k+1) - lgamma(a+1) - lgamma(b+1) - lgamma(c+1) + (a+b)*log(.25) + c*log(.5)
        lpD.setdefault(a-b, []).append(x)
Ds = sorted(lpD); lpN = [lse(lpD[t]) for t in Ds]
lcN = [lb(k, j) + j*log(3) for j in range(k+1)]   # register projection: Paulis on <= w sites
cN = windows(lpN, lcN)
nat_up = [(2*q/(3*L2), -e/(3*L2)) for q, e in cN]      # ring <= |pi|^2
nat_lo = [(q/(2*L2), -e/(3*L2)) for q, e in cN]        # ring >= |pi|^{3/2} (optimistic)
# (c) S'_3 bound, d=2: Pr[L1 = k-j] = f(j) * (k-2j+1)/2^k, f(j) = C(k,j)(k-2j+1)/(k-j+1)
lpL = {}
for j in range(k//2+1):
    lpL[k-j] = lb(k, j) + 2*log(k-2*j+1) - log(k-j+1) - k*log(2)
Lvals = sorted(lpL)
import numpy as np
P = np.array([math.exp(lpL[l]) for l in Lvals]); Lv = np.array(Lvals)
cdf = np.cumsum(P)  # Pr[L <= Lv[i]]
def pr_gap(u):  # Pr[L1 > max(L2,L3) + u] = sum_l P(l) * Pr[L <= l-u-1]^2  (computed in logs crudely)
    tot = -1e301
    for i, l in enumerate(Lvals):
        m = l - u - 1
        idx = np.searchsorted(Lv, m, side='right') - 1
        if idx < 0: continue
        tot = lse([tot, lpL[l] + 2*math.log(max(cdf[idx], 1e-300))])
    return tot
def u_of_q(lq):  # max j with C(k,j)(k-2j+1)/(k-j+1) <= q  (two-row minimiser, as in S')
    j = 0
    while j+1 <= k//2 and lb(k, j+1) + log(k-2*(j+1)+1) - log(k-(j+1)+1) <= lq: j += 1
    return j
def interp(curve, v):
    best = 0.0
    for vv, du in curve:
        if vv <= v + 1e-12: best = max(best, du)
    return best
print(f'k={k}   columns: v | composite(Majorana) 1-u | native upper-ring | native optimistic-ring | S3-bound (max 1-u for symmetric)')
for v in [0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3]:
    a = interp(compo, v); b = interp(nat_up, v); c = interp(nat_lo, v)
    lq = 3*v*L2; u = u_of_q(lq); bnd = -pr_gap(u)/(3*L2)
    print(f'v={v:5.2f} | comp {a:.5f} ({a/v**2:.3f} v^2) | nat {b:.5f} ({b/v**2:.3f}) | natopt {c:.5f} ({c/v**2:.3f}) | S3 <= {bnd:.5f} ({bnd/v**2:.3f} v^2, u={u})')
