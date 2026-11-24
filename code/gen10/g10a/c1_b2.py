# Independent check of Theorem B2's combinatorial core and constants.
import itertools, random, math
from math import comb, sqrt
# (a) E|B-B'| = N*C(2N,N)/4^N ; ratio to sqrt(N)
mn = min((N*comb(2*N,N)/4**N)/sqrt(N) for N in range(1,400))
print("min_N E|B-B'|/sqrt(N) =", mn, " (N=1,2,3:", [N*comb(2*N,N)/4**N for N in (1,2,3)], ")")
print("c_3 =", 0.189/sqrt(2*comb(3,2)))
# (b) core: for any family (partition_t, alpha_t), t in X1:
#   max_omega sum_t def^omega_t >= 0.189 * sum_v sqrt(N_v)
random.seed(1)
def trial(n, m, T, g, mode):
    # random 3-CNF-like clause list with degree <= 3
    occ = {v: [] for v in range(n)}
    clauses = []
    pool = [v for v in range(n) for _ in range(3)]
    random.shuffle(pool)
    while len(pool) >= 3 and len(clauses) < m:
        e = set()
        tries = 0
        while len(e) < 3 and pool and tries < 50:
            v = pool.pop(); tries += 1
            if v in e: pool.insert(0, v)
            else: e.add(v)
        if len(e) == 3:
            for v in e: occ[v].append(len(clauses))
            clauses.append(tuple(e))
    designs = []
    for t in range(T):
        lab = [random.randrange(g) for _ in range(n)]
        if mode == 0:   alpha = [random.randrange(g) for _ in clauses]
        elif mode == 1: alpha = [lab[random.choice(e)] for e in clauses]
        else:           # adversarial-ish: try to keep variables unsplit
            alpha = [lab[min(e)] for e in clauses]
        designs.append((lab, alpha))
    # N_v
    S = 0.0
    for v in range(n):
        best = 0
        for o1, o2 in itertools.combinations(occ[v], 2):
            best = max(best, sum(1 for lab, al in designs if al[o1] != al[o2]))
        S += sqrt(best)
    # brute-force max over omega of sum_t def
    choices = [occ[v] if occ[v] else [None] for v in range(n)]
    best = -1; worst = 10**9
    for om in itertools.product(*choices):
        tot = 0
        for lab, al in designs:
            a = [0]*g; X = [0]*g
            for v in range(n):
                a[lab[v]] += 1
                X[al[om[v]] if om[v] is not None else lab[v]] += 1
            tot += sum(max(a[h]-X[h], 0) for h in range(g))
        best = max(best, tot); worst = min(worst, tot)
    return best, 0.189*S, worst
minratio = 1e9
for it in range(60):
    n = random.choice([6,7,8]); T = random.choice([1,2,3,5,8]); g = random.choice([2,3,4])
    b, lb, w = trial(n, n, T, g, it % 3)
    if lb > 0: minratio = min(minratio, b/lb)
    assert b >= lb - 1e-9, (b, lb)
print("B2 core: 60 brute-force instances pass; min (max_omega sum def)/(0.189 sum sqrt N_v) =", round(minratio,3))
