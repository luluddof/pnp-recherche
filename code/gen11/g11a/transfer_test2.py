# Clustered instances where small genuine separators exist and the lemma's bounds are non-trivial (small L).
import random, itertools
from fractions import Fraction as Fr
import importlib.util, sys
src = open('transfer_test.py').read().replace("\nmain()\n", "\n")
exec(src)
def expand(n, clauses, rnd, keep=.5):
    m = len(clauses); d = [0]*n
    for e in clauses:
        for v in e: d[v] += 1
    pi = []; fib = {v: [] for v in range(n)}; sc = []
    for e in clauses:
        cs = []
        for v in e:
            x = len(pi); pi.append(v); fib[v].append(x); cs.append(x)
        sc.append(tuple(cs))
    for v in range(n):
        f = fib[v][:]; rnd.shuffle(f)
        if len(f) >= 2:
            for i in range(len(f)): sc.append((f[i], f[(i+1) % len(f)]))
    n0 = 0
    for v in range(n):
        if d[v] == 0 and rnd.random() < keep:
            x = len(pi); pi.append(v); fib[v].append(x); n0 += 1
    ns = len(pi); adj = [set() for _ in range(ns)]
    for c in sc:
        for a in c:
            for b in c:
                if a != b: adj[a].add(b)
    occ = [0]*ns
    for c in sc:
        for x in set(c): occ[x] += 1
    assert max(occ) <= 3
    return d, pi, fib, adj, ns, n0
rnd = random.Random(7)
st = dict(sep=0, checked=0, nontrivialA=0, nontrivialB=0, minslackA=10**9, minslackB=10**9, inst=0)
for it in range(150):
    nb = rnd.randint(2, 3); per = rnd.randint(4, 6); mper = rnd.randint(3, 6)
    n = nb*per + 2; clauses = []
    for b in range(nb):
        for _ in range(mper):
            clauses.append(tuple(x + b*per for x in rnd.sample(range(per), 3)))
    for _ in range(rnd.randint(0, 2)):
        b1, b2 = rnd.sample(range(nb), 2)
        clauses.append((b1*per + rnd.randrange(per), b2*per + rnd.randrange(per), n-1))
    d, pi, fib, adj, ns, n0 = expand(n, clauses, rnd)
    st['inst'] += 1
    for k in range(0, min(ns//9, 2)+1):
        for Ys in itertools.combinations(range(ns), k):
            check_one(n, len(clauses), clauses, d, pi, fib, adj, ns, n0, Ys, st)
print("clustered:", st)
