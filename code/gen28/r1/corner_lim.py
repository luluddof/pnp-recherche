# k -> infinity limit: eta = H(tau) + tau log2 3 ; delta = max over window [a, a+tau] of min(D(a||1/4), D(a+tau||1/4)) (bits)
from math import log2, log
from scipy.optimize import brentq
H = lambda t: -t*log2(t) - (1-t)*log2(1-t)
def D(a, p=0.25):
    if a <= 0: return log2(1/(1-p))
    return a*log2(a/p) + (1-a)*log2((1-a)/(1-p))
for tau in [1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-8, 1e-12]:
    g = lambda a: D(a) - D(a+tau)
    a = brentq(g, max(1e-15, 0.25 - tau), 0.25)
    dl = D(a); eta = H(tau) + tau*log2(3)
    print(f"tau={tau:.0e} eta={eta:.3e} delta={dl:.3e} delta/eta^2={dl/eta**2:.4f} delta*log2(1/eta)^2/eta^2={dl*log2(1/eta)**2/eta**2:.4f}")
print("predicted limit 2/(3 ln2) =", 2/(3*log(2)))
