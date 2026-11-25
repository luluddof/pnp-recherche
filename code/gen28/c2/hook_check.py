# check f^mu >= C(k,j)/sqrt(e) for j = k - mu_1 <= k/4 (all partitions, k <= 26)
from sympy.utilities.iterables import partitions
from math import comb, factorial, sqrt, e
worst = 9
for k in range(2, 27):
    for p in partitions(k):
        mu = sorted([a for a,c in p.items() for _ in range(c)], reverse=True)
        j = k - mu[0]
        if 4*j > k: continue
        cols = [sum(1 for r in mu if r > c) for c in range(mu[0])]
        h = 1
        for i,r in enumerate(mu):
            for c in range(r): h *= (r-c-1) + (cols[c]-i-1) + 1
        f = factorial(k)//h
        worst = min(worst, f*sqrt(e)/comb(k,j))
print("min over k<=26, j<=k/4 of f^mu*sqrt(e)/C(k,j) =", worst)
