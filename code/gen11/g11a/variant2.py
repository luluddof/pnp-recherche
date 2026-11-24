import random, itertools
from fractions import Fraction as Fr
exec(open('variant.py').read().split("rnd = random.Random(11)")[0])
rnd = random.Random(5); checked = 0; minA = minB = 10**9; nontriv = 0
for it in range(200):
    nb = rnd.randint(2, 3); per = rnd.randint(4, 7); n = nb*per + 1; clauses = []
    for b in range(nb):
        for _ in range(rnd.randint(2, 7)): clauses.append(tuple(x + b*per for x in rnd.sample(range(per), 3)))
    for _ in range(rnd.randint(0, 2)):
        b1, b2 = rnd.sample(range(nb), 2); clauses.append((b1*per + rnd.randrange(per), b2*per + rnd.randrange(per), n-1))
    d, pi, fib, adj, ns = expand_var(n, clauses, rnd)
    for k in range(0, min(ns//9, 2) + 1):
        for Ys in itertools.combinations(range(ns), k):
            comps = components(ns, adj, Ys)
            if any(3*len(c) > 2*ns for c in comps): continue
            Y = set(pi[x] for x in Ys); L = sum(d[v] for v in Y); kk = len(comps)
            if kk > 10: continue
            for mask in range(1 << kk):
                sz = sum(len(comps[i]) for i in range(kk) if mask >> i & 1)
                if not (3*sz >= ns and 3*sz <= 2*ns): continue
                As = set(x for i in range(kk) if mask >> i & 1 for x in comps[i])
                A = set(v for v in range(n) if v not in Y and fib[v] and all(x in As for x in fib[v]))
                B = set(range(n)) - A - Y
                mA = sum(all(v in A for v in e) for e in clauses); mB = sum(all(v in B for v in e) for e in clauses)
                for e in clauses: assert all(v in A for v in e) or all(v in B for v in e) or any(v in Y for v in e)
                assert mA >= Fr(ns, 9) - L and mB >= Fr(2*ns, 27) - L
                nontriv += (Fr(2*ns, 27) - L > 0)
                minA = min(minA, mA - (Fr(ns, 9) - L)); minB = min(minB, mB - (Fr(2*ns, 27) - L)); checked += 1
print("variant (copy only d_v>3) checked", checked, "nontrivial B-bound", nontriv, "min slack", minA, minB)
