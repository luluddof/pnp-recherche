# G10c: build the real 3-CNF, count its models exactly with an independent DPLL counter (count.cpp), and compare with
#   sum over admissible (Y,X) of prod_{e in X} w_e * p^{cc(Y,X)}   (brute force from the definition),
# then check (#Psi/p) mod p against the independently computed (weighted) number of trees / cycles / paths.
import itertools, random, subprocess, sys
sys.path.insert(0, '.')
from identity_all import cc_of, tau, count_cycles, count_paths
from encode import build, dimacs
rnd = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
def nmodels(F):
    return int(subprocess.run(['./count'], input=dimacs(F), capture_output=True, text=True, check=True).stdout)
def direct(n, E, p, w, kind, l, K):
    tot = 0; unit = 0
    Kmask = sum(1 << v for v in K)
    for Y in range(1 << n):
        if Y & Kmask != Kmask: continue
        ny = bin(Y).count('1')
        idx = [j for j, e in enumerate(E) if Y >> e[0] & 1 and Y >> e[1] & 1 and w[j] != 0]
        if kind == 'ST':
            if ny != l + 1: continue
            r = l
        elif kind == 'CYC':
            if ny != l: continue
            r = None
        else:
            if ny != l: continue
            r = l - 1
        for rr in ([r] if r is not None else range(len(idx) + 1)):
            for Xi in itertools.combinations(idx, rr):
                deg = [0] * n
                for j in Xi: deg[E[j][0]] += 1; deg[E[j][1]] += 1
                if kind == 'CYC' and not all(deg[v] == 2 for v in range(n) if Y >> v & 1): continue
                if kind == 'PATH' and max(deg) > 2: continue
                c = cc_of(Y, [E[j] for j in Xi], n); pr = 1
                for j in Xi: pr *= w[j]
                tot += pr * p ** c; unit += 1
    return tot, unit
tests = 0; maxratio = 0; maxw = 0
for trial in range(int(sys.argv[1])):
    n = rnd.randint(3, 6); m = rnd.randint(n - 1, min(n * (n - 1) // 2, 8))
    allE = list(itertools.combinations(range(n), 2)); rnd.shuffle(allE); E = sorted(allE[:m])
    p = rnd.choice([2, 3, 3, 5, 5, 6, 7])
    mode = rnd.choice(['ones', 'rand', 'randp'])
    w = [1] * m if mode == 'ones' else [rnd.randrange(p) for _ in E] if mode == 'rand' else [rnd.randint(1, p) for _ in E]
    kind = rnd.choice(['ST', 'ST', 'CYC', 'PATH'])
    K = rnd.sample(range(n), rnd.randint(1, 3)) if kind == 'ST' else []
    l = rnd.randint(1, n - 1) if kind == 'ST' else rnd.randint(3, n) if kind == 'CYC' else rnd.randint(2, n)
    F, M0, b = build(n, E, p, w, kind, l, K)
    assert all(len(c) <= 3 for c in F.cl); maxw = max(maxw, max(len(c) for c in F.cl))
    maxratio = max(maxratio, len(F.cl) / (b * (n + m)))
    got = nmodels(F); want, npairs = direct(n, E, p, w, kind, l, K)
    assert got == want, (kind, n, E, p, w, l, K, got, want)
    # skeleton alone (p-part removed): #models == number of admissible (x,y)  [unique extension, property (ii)]
    # check mod-p extraction against an independent count
    adj = [set() for _ in range(n)]
    for j, (u, v) in enumerate(E):
        if w[j] != 0: adj[u].add(v); adj[v].add(u)
    if mode == 'ones':
        if kind == 'ST':
            Km = sum(1 << v for v in K)
            ref = sum(tau(Yt, E) for Yt in itertools.combinations(range(n), l + 1) if sum(1 << v for v in Yt) & Km == Km)
        elif kind == 'CYC': ref = count_cycles(n, adj, l)
        else: ref = count_paths(n, adj, l)
        assert got % p == 0 and (got // p) % p == ref % p, (kind, got, ref)
    elif kind == 'ST' and mode == 'rand':
        Km = sum(1 << v for v in K)
        ref = sum(tau(Yt, E, wt=w) for Yt in itertools.combinations(range(n), l + 1) if sum(1 << v for v in Yt) & Km == Km)
        assert got % p == 0 and (got // p) % p == ref % p, (kind, got, ref)
    tests += 1
print(f"CNF check: {tests} formulas, #Psi == direct weighted sum in all; max clause width {maxw}; max clauses/(b(n+m)) = {maxratio:.1f}")
