# Theorem S' quantitative curve at d=2 (qubit sites), N = 2^k:
#   S_k-invariant K, osr(K) <= q  =>  u <= u(q) := max{k-mu_1 : mu |- k, l(mu) <= 4, f^mu <= q}
#   p >= N^2 * Pr[m2 - m1 > 2u]   (m = 2nd row of SW shape of a uniform binary word; m1, m2 iid)
# Compare with the best S_k-symmetric construction known (T3 window family: K = f(M), deg f = w).
import numpy as np, sys
from math import lgamma, log, comb, exp
LN2 = log(2)

def logf4(mu):  # mu: 4 parts (with zeros), log of number of SYT, via f = k! prod_{i<j}(l_i-l_j) / prod l_i!
    r = len(mu); l = [mu[i] + (r-1-i) for i in range(r)]; k = sum(mu)
    s = lgamma(k+1)
    for i in range(r):
        for j in range(i+1, r): s += log(l[i]-l[j])
    for i in range(r): s -= lgamma(l[i]+1)
    return s

def g_table(k, jmax):
    g = []
    for j in range(jmax+1):
        best = None; m1 = k - j
        for a in range(min(j, m1), -1, -1):
            for b in range(min(j-a, a), -1, -1):
                c = j - a - b
                if c < 0 or c > b: continue
                v = logf4([m1, a, b, c])
                best = v if best is None else min(best, v)
        g.append(best)
    return g

def sw_logpr(k):
    out = []
    for m in range(k//2+1):
        f = comb(k,m) - (comb(k,m-1) if m>0 else 0)
        out.append(log(f) + log(k-2*m+1) - k*LN2)
    return np.array(out)

def tail_bound(k, u, lp):
    pr = np.exp(lp); M = len(pr)
    cdf_tail = np.concatenate([np.cumsum(pr[::-1])[::-1], [0.0]])  # Pr[m >= t]
    tot = 0.0
    for a in range(M):
        t = a + 2*u + 1
        if t < M: tot += pr[a]*cdf_tail[t]
    return tot

def family(k):
    # T3 window family: list of (log q, log p) over w
    res = []
    for w in range(0, k):
        q = sum(comb(k,j)*3**j for j in range(w+1))
        best = None
        for lo in range(0, k-w+1):
            p = sum(comb(k,m)*3**(k-m) for m in range(k+1) if not (lo <= m <= lo+w))
            best = p if best is None else min(best, p)
        if best == 0: break
        res.append((log(q), log(best)))
    return res

for k in [int(a) for a in sys.argv[1:]]:
    lnN = k*LN2
    g = g_table(k, k//4)
    lp = sw_logpr(k)
    print(f"== k={k} (N=2^{k})")
    print("  S' bound:   eta     u   delta_max   delta/eta^2   delta*log2(1/eta)^2/eta^2")
    for eta in [0.02, 0.05, 0.1, 0.2, 0.3, 0.5]:
        lq = eta*lnN
        u = max([j for j in range(len(g)) if g[j] <= lq + 1e-9], default=0)
        if u >= k//4 - 1: print(f"   eta={eta}: u hits k/4 cap"); continue
        b = tail_bound(k, u, lp)
        if b <= 0: print(f"   eta={eta}: bound 0"); continue
        delta = -log(b)/lnN
        print(f"   {eta:5.2f}  {u:4d}   {delta:.5f}    {delta/eta**2:.3f}      {delta*np.log2(1/eta)**2/eta**2:.3f}")
    if k <= 400:
        fam = family(k)
        print("  family (T3 window): eta, delta, delta/eta^2")
        for lq, lpp in fam[1:8]:
            eta = lq/lnN; delta = 2 - lpp/lnN
            print(f"   {eta:.3f}  {delta:.5f}  {delta/eta**2:.3f}")
    sys.stdout.flush()
