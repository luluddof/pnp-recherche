# is min_{mu |- k, mu_1 = k-j} f^mu attained at (k-j, j) for j <= k/2 ?  (all rows allowed)
from sympy.utilities.iterables import partitions
from math import factorial
bad = 0
for k in range(2, 31):
    best = {}
    for p in partitions(k):
        mu = sorted([a for a,c in p.items() for _ in range(c)], reverse=True)
        j = k - mu[0]
        if 2*j > k: continue
        cols = [sum(1 for r in mu if r > c) for c in range(mu[0])]
        h = 1
        for i,r in enumerate(mu):
            for c in range(r): h *= (r-c-1) + (cols[c]-i-1) + 1
        f = factorial(k)//h
        if j not in best or f < best[j][0]: best[j] = (f, tuple(mu))
    for j,(f,mu) in best.items():
        if j>0 and len(mu) != 2: bad += 1; print("non-two-row minimiser", k, j, mu)
print("violations:", bad)
