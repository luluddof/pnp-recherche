# Variant "copy only variables with d_v > 3" (fibre size 1 if 1<=d_v<=3): generalized transfer bound with n* in place of 3m
import random, itertools, math
from fractions import Fraction as Fr
src = open('transfer_test.py').read().replace("\nmain()\n", "\n"); exec(src)
def expand_var(n, clauses, rnd):
    m = len(clauses); d = [0]*n
    for e in clauses:
        for v in e: d[v] += 1
    pi = []; fib = {v: [] for v in range(n)}; sc = []
    single = {}
    for e in clauses:
        cs = []
        for v in e:
            if d[v] <= 3:
                if v not in single:
                    single[v] = len(pi); pi.append(v); fib[v].append(single[v])
                cs.append(single[v])
            else:
                x = len(pi); pi.append(v); fib[v].append(x); cs.append(x)
        sc.append(tuple(cs))
    for v in range(n):
        f = fib[v][:]; rnd.shuffle(f)
        if len(f) >= 2:
            for i in range(len(f)): sc.append((f[i], f[(i+1) % len(f)]))
    ns = len(pi); adj = [set() for _ in range(ns)]; occ = [0]*ns
    for c in sc:
        for x in set(c): occ[x] += 1
        for a in c:
            for b in c:
                if a != b: adj[a].add(b)
    assert max(occ) <= 3
    return d, pi, fib, adj, ns
rnd = random.Random(11); checked = 0; minA = minB = 10**9
for it in range(300):
    n = rnd.randint(4, 9); m = rnd.randint(6, 14)
    clauses = [tuple(rnd.sample(range(n), 3)) for _ in range(m)]
    d, pi, fib, adj, ns = expand_var(n, clauses, rnd)
    for k in range(0, ns//9 + 1):
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
                # generalized: mA >= n*/9 - L, mB >= 2n*/27 - L  (n0 = 0 here)
                assert mA >= Fr(ns, 9) - L and mB >= Fr(2*ns, 27) - L
                minA = min(minA, mA - (Fr(ns, 9) - L)); minB = min(minB, mB - (Fr(2*ns, 27) - L)); checked += 1
print("variant checked", checked, "min slack", minA, minB)
H = lambda x: -x*math.log(x)-(1-x)*math.log(1-x)
kap = lambda b, e: 3*(1-e)*H(b)-math.log(2)-H(e)
# xi for the variant at density Delta (Poisson(3 Delta) degrees): n*/(3m) -> 1 - [P(2)+2P(3)]/(3 Delta)
for D in [2, 3, 4.26]:
    lam = 3*D; P = lambda k: math.exp(-lam)*lam**k/math.factorial(k)
    xi = 1 - (P(2) + 2*P(3))/lam
    print("Delta", D, "xi", xi, "2xi/9 - 1/20 - 1/360 =", 2*xi/9 - 1/20 - 1/360, "(needs >= beta)")
print("kappa(0.16,0.05)", kap(.16, .05), "threshold", math.log(2)/kap(.16, .05))
print("note-4 numbers: xi=1-1/Delta; Delta=2 ->", 1-1/2, "; needs xi>2/3 i.e. Delta>3")
