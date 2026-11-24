# Independent exhaustive test of the Transfer lemma (G11a).
import random, itertools, sys
from fractions import Fraction as Fr

def build(n, m, rnd, keep_prob=0.5, variant_tovey=True):
    clauses = [tuple(rnd.sample(range(n), 3)) for _ in range(m)]
    d = [0]*n
    for e in clauses:
        for v in e: d[v] += 1
    pi = []            # copy -> original variable
    fib = {v: [] for v in range(n)}
    star_clauses = []  # list of variable sets (copies)
    for e in clauses:
        cs = []
        for v in e:
            x = len(pi); pi.append(v); fib[v].append(x); cs.append(x)
        star_clauses.append(tuple(cs))
    for v in range(n):
        f = fib[v][:]; rnd.shuffle(f)
        s = len(f)
        if s >= 2:
            for i in range(s):
                a, b = f[i], f[(i+1) % s]
                star_clauses.append((a, b))   # (not a or b)
    n0 = 0
    for v in range(n):
        if d[v] == 0 and rnd.random() < keep_prob:
            x = len(pi); pi.append(v); fib[v].append(x); n0 += 1
    ns = len(pi)
    assert ns == 3*m + n0
    # occurrences per copy <= 3
    occ = [0]*ns
    for c in star_clauses:
        for x in set(c): occ[x] += 1
    assert max(occ) <= 3, max(occ)
    adj = [set() for _ in range(ns)]
    for c in star_clauses:
        for a in c:
            for b in c:
                if a != b: adj[a].add(b)
    return clauses, d, pi, fib, adj, ns, n0

def components(ns, adj, Ys):
    seen = set(Ys); comps = []
    for s in range(ns):
        if s in seen: continue
        comp = [s]; seen.add(s); st = [s]
        while st:
            x = st.pop()
            for y in adj[x]:
                if y not in seen: seen.add(y); comp.append(y); st.append(y)
        comps.append(comp)
    return comps

def check_one(n, m, clauses, d, pi, fib, adj, ns, n0, Ys, stats):
    comps = components(ns, adj, Ys)
    if any(3*len(c) > 2*ns for c in comps): return  # not a 2/3-separator
    if 9*len(Ys) > ns: return
    stats['sep'] += 1
    Y = set(pi[x] for x in Ys)
    assert len(Y) <= len(Ys)
    L = sum(d[v] for v in Y)
    # all unions A* of components with n*/3 <= |A*| <= 2n*/3 (the lemma only needs one; test all, <= 2^k)
    k = len(comps)
    cand = []
    if k <= 12:
        for mask in range(1 << k):
            sz = sum(len(comps[i]) for i in range(k) if mask >> i & 1)
            if 3*sz >= ns and 3*sz <= 2*ns: cand.append(mask)
    else:
        # greedy as in the proof
        comps.sort(key=len, reverse=True)
        if 3*len(comps[0]) >= ns: cand = [1]
        else:
            mask = 0; sz = 0; i = 0
            while 3*sz < ns: mask |= 1 << i; sz += len(comps[i]); i += 1
            assert 3*sz <= 2*ns
            cand = [mask]
    assert cand, "no A* exists"
    for mask in cand:
        As = set(x for i in range(k) if mask >> i & 1 for x in comps[i])
        Bs = set(range(ns)) - As - set(Ys)
        assert 9*len(Bs) >= 2*ns
        A = set(); B = set()
        for v in range(n):
            if v in Y: continue
            f = fib[v]
            if f and all(x in As for x in f): A.add(v)
            else:
                assert all(x in Bs for x in f), "fibre split (E1 violated?)"
                B.add(v)
        mA = mB = cross = meetY = 0
        for e in clauses:
            inA = all(v in A for v in e); inB = all(v in B for v in e)
            hasY = any(v in Y for v in e)
            assert inA or inB or hasY, "claim (1) fails"
            mA += inA; mB += inB; meetY += hasY
            if any(v in A for v in e) and any(v not in A for v in e): cross += 1
        assert cross <= L, "crossing > L"
        assert meetY <= L
        lhsA = Fr(m, 3) - L - Fr(n0, 3)
        lhsB = Fr(2*m, 9) - L - Fr(n0, 3)
        assert mA >= lhsA, ("claim (2) A fails", mA, lhsA)
        assert mB >= lhsB, ("claim (2) B fails", mB, lhsB)
        stats['checked'] += 1
        if lhsA > 0: stats['nontrivialA'] += 1
        if lhsB > 0: stats['nontrivialB'] += 1
        stats['minslackA'] = min(stats['minslackA'], mA - lhsA)
        stats['minslackB'] = min(stats['minslackB'], mB - lhsB)

def main():
    stats = dict(sep=0, checked=0, nontrivialA=0, nontrivialB=0, minslackA=10**9, minslackB=10**9, inst=0)
    rnd = random.Random(20261124)
    # exhaustive over all Y* with |Y*| <= n*/9
    for it in range(400):
        n = rnd.randint(4, 14); m = rnd.randint(3, 12)
        inst = build(n, m, rnd, keep_prob=rnd.choice([0, .5, 1]))
        clauses, d, pi, fib, adj, ns, n0 = inst
        stats['inst'] += 1
        kmax = ns // 9
        for k in range(0, kmax+1):
            for Ys in itertools.combinations(range(ns), k):
                check_one(n, m, clauses, d, pi, fib, adj, ns, n0, Ys, stats)
    print("exhaustive:", stats)
    # structured instances: disjoint blocks so that small separators exist and bounds bite (L small)
    stats2 = dict(sep=0, checked=0, nontrivialA=0, nontrivialB=0, minslackA=10**9, minslackB=10**9, inst=0)
    for it in range(300):
        # two or three clusters joined through few shared variables
        nb = rnd.randint(2, 4); per = rnd.randint(4, 7); mper = rnd.randint(4, 9)
        n = nb*per + 2
        clauses = []
        for b in range(nb):
            for _ in range(mper):
                clauses.append(tuple(x + b*per for x in rnd.sample(range(per), 3)))
        for _ in range(rnd.randint(0, 2)):
            # bridge clause through the hub variable n-1
            b1, b2 = rnd.sample(range(nb), 2)
            clauses.append((b1*per + rnd.randrange(per), b2*per + rnd.randrange(per), n-1))
        m = len(clauses)
        # reuse build's logic by monkeypatching the clause list
        class R2(random.Random): pass
        d = [0]*n
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
            if d[v] == 0 and rnd.random() < .5:
                x = len(pi); pi.append(v); fib[v].append(x); n0 += 1
        ns = len(pi)
        adj = [set() for _ in range(ns)]
        for c in sc:
            for a in c:
                for b in c:
                    if a != b: adj[a].add(b)
        stats2['inst'] += 1
        kmax = min(ns // 9, 3)
        for k in range(0, kmax+1):
            for Ys in itertools.combinations(range(ns), k):
                check_one(n, m, clauses, d, pi, fib, adj, ns, n0, Ys, stats2)
    print("clustered :", stats2)
main()
