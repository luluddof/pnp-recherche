# generalised Kronecker coefficients g(l1,l2,l3,l4) = multiplicity of S_l1 (x) .. (x) S_l4 in S^d(V^{(x)4})
from functools import lru_cache
from math import factorial
from collections import Counter
from monos import parts
def partitions(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield (); return
    for k in range(min(n, mx), 0, -1):
        for r in partitions(n-k, k): yield (k,) + r
@lru_cache(None)
def chi(lam, mu):   # Murnaghan-Nakayama, lam partition (tuple, no zeros), mu cycle type (tuple)
    if not mu: return 1 if sum(lam) == 0 else 0
    k = mu[0]; rest = mu[1:]
    # beta numbers
    n = len(lam); beta = [lam[i] + n - 1 - i for i in range(n)]; s = set(beta); tot = 0
    for b in beta:
        if b - k >= 0 and (b - k) not in s:
            nb = sorted([x for x in beta if x != b] + [b - k], reverse=True)
            sign = (-1) ** sum(1 for x in beta if b - k < x < b)
            lam2 = tuple(x for x in (nb[i] - (n - 1 - i) for i in range(n)) if x > 0)
            tot += sign * chi(lam2, rest)
    return tot
def classsize(mu):
    d = sum(mu); c = Counter(mu); den = 1
    for k, m in c.items(): den *= k**m * factorial(m)
    return factorial(d) // den
def g(*lams):
    lams = [tuple(x for x in l if x > 0) for l in lams]; d = sum(lams[0]); tot = 0
    for mu in partitions(d):
        pr = classsize(mu)
        for l in lams: pr *= chi(l, mu)
        tot += pr
    assert tot % factorial(d) == 0
    return tot // factorial(d)
if __name__ == "__main__":
    for d in range(3, 13):
        P = [m for m in parts(d) if m[2] > 0]
        row = []
        for mu in P:
            for nu in P:
                if nu <= mu: row.append((mu, nu, g(mu, nu, nu, mu)))
        print(d, row, flush=True)
