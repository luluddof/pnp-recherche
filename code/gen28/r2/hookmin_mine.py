from math import factorial, comb
def parts(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield (); return
    for a in range(min(n, mx), 0, -1):
        for r in parts(n-a, a): yield (a,) + r
def f(mu):
    k = sum(mu); cols = [sum(1 for r in mu if r > c) for c in range(mu[0])]; h = 1
    for i, r in enumerate(mu):
        for c in range(r): h *= (r-c-1) + (cols[c]-i-1) + 1
    return factorial(k)//h
bad = 0
for k in range(2, 29):
    best = {}
    for mu in parts(k):
        j = k - mu[0]
        if j > k//2: continue
        v = f(mu)
        if j not in best or v < best[j][0]: best[j] = (v, mu)
    for j, (v, mu) in best.items():
        two = f((k-j, j)) if j else 1
        if v < two: bad += 1; print("counterex", k, j, mu, v, two)
    # hook lemma check
    for mu in parts(k):
        j = k-mu[0]
        if 4*j <= k and f(mu)*1.6487212707 < comb(k, j): bad += 1; print("hook lemma fails", mu)
print("violations", bad)
