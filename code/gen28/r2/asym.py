# d=2 exact finite-k S' bound (C2 threshold 2u vs improved threshold u) and asymptotic formulas.
import numpy as np
from scipy.special import gammaln
from scipy.optimize import brentq
def lc(n, r): return gammaln(n+1)-gammaln(r+1)-gammaln(n-r+1)
def run(k, eta):
    ms = np.arange(k//2+1)
    logf = np.array([lc(k,m) + np.log((k-2*m+1)/(k-m+1)) for m in ms])  # f^{(k-m,m)}
    logpr = logf + np.log(k-2*ms+1) - k*np.log(2)
    pr = np.exp(logpr - logpr.max()); pr /= pr.sum()  # tiny renorm issue only
    logq = eta*k*np.log(2)
    u = max(m for m in ms if logf[m] <= logq + 1e-12)   # two-row minimiser (C2 hook_min, j<=k/2)
    L = k - ms
    out = []
    for thr in [2*u, u]:
        # Pr[L1-L2>thr] = sum_{a,b} pr[a]pr[b] [L[a]-L[b]>thr]
        cdf = np.cumsum(pr[::-1])[::-1]  # Pr[m >= b]
        P = 0.0
        for a in ms:
            # L[a]-L[b] = b - a > thr -> b > a+thr
            b0 = a + thr + 1
            if b0 <= k//2: P += pr[a]*cdf[b0]
        out.append(-np.log2(P)/k if P > 0 else np.inf)
    return u, out
def asymp(eta, c):
    return c*eta**2/np.log(1/eta)**2
for k in [400, 1600, 6400]:
    for eta in [0.02, 0.05, 0.1, 0.2, 0.3]:
        u, (d2u, du) = run(k, eta)
        print(f"k={k} eta={eta}: u={u} delta_bound(2u)={d2u:.4f} delta_bound(u)={du:.4f}  formula 5.545:{asymp(eta,5.545):.4f} 1.386:{asymp(eta,1.386):.4f}")
# exact k->oo limit: x solves H_e(x)=eta ln2 ; rate D(1/2+t||1/2)/ln2, t = 2x or x
H = lambda x: -x*np.log(x)-(1-x)*np.log(1-x)
D = lambda t: (0.5+t)*np.log(1+2*t)+(0.5-t)*np.log(1-2*t)
for eta in [1e-1, 1e-2, 1e-3, 1e-4, 1e-6]:
    x = brentq(lambda x: H(x)-eta*np.log(2), 1e-15, 0.5)
    lim2, lim1 = D(min(2*x,0.4999))/np.log(2), D(x)/np.log(2)
    print(f"eta={eta:g}: limit delta(2u)={lim2:.3e} ratio to 5.545 formula {lim2/asymp(eta,5.545):.3f}; delta(u)={lim1:.3e} ratio to 1.386 formula {lim1/asymp(eta,1.386):.3f}")
