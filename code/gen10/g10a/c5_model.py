import sys, random, math, itertools
sys.path.insert(0, '.')
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from c4_thmC import wtw, rank_mod, P
# (A) group-atomic charge W_t vs row-level variable elimination.
# k measured groups, r rows each; clauses assigned to group i read group j only through ONE functional (row i of group j).
# Factor of group i:  N_i( all rows of i ; row (j,i) for j != i ).
for k, r in ((3, 4), (3, 5), (4, 4)):
    V = [(i, j) for i in range(k) for j in range(r)]; ix = {v: n for n, v in enumerate(V)}
    adj = {n: set() for n in range(len(V))}
    for i in range(k):
        scope = [(i, j) for j in range(r)] + [(j, i) for j in range(k) if j != i]
        for x in scope:
            for y in scope:
                if x != y: adj[ix[x]].add(ix[y])
    row_w = wtw(len(V), adj, [1]*len(V))
    gadj = {i: set(range(k)) - {i} for i in range(k)}
    W = wtw(k, gadj, [r]*k)
    print(f"k={k} groups x r={r} rows: G7b charge W_t = {W};  row-level elimination: max bag = {row_w} rows  (formula max(k(k-1), r+k-1) = {max(k*(k-1), r+k-1)})")
# (B) Lemma E numerically: E(a-r-X)_+ >= (c0/2) min(V, sqrt V) - r for Poisson-binomial X with mean <= a
rng = random.Random(5); c0 = 0.216; worst = 1e9
for it in range(4000):
    m = rng.randint(1, 40)
    p = [rng.choice([1/3, 2/3, 1/2, 1.0, rng.random()]) for _ in range(m)]
    dist = [1.0]
    for q in p:
        nd = [0.0]*(len(dist)+1)
        for i, x in enumerate(dist): nd[i] += x*(1-q); nd[i+1] += x*q
        dist = nd
    mu = sum(p); Vr = sum(q*(1-q) for q in p)
    a = math.ceil(mu) + rng.randint(0, 3); r = rng.randint(0, 3)
    lhs = sum(pr*max(a-r-x, 0) for x, pr in enumerate(dist))
    rhs = c0/2*min(Vr, math.sqrt(Vr)) - r
    worst = min(worst, lhs - rhs); assert lhs >= rhs - 1e-12
    e_abs = sum(pr*abs(a-x) for x, pr in enumerate(dist))
    assert e_abs >= c0*min(Vr, math.sqrt(Vr)) - 1e-12
print("Lemma E / anticoncentration: 4000 Poisson-binomial cases pass, min slack", round(worst, 4))
# (C) clean generic test: complete interaction, generic rows in U_o => nullity == N^omega exactly
rng = random.Random(11); eq = tot = 0
for it in range(200):
    n = rng.choice([8, 10, 12, 14]); g = rng.choice([2, 3, 4])
    lab = [rng.randrange(g) for _ in range(n)]
    S = [[v for v in range(n) if lab[v] == h] for h in range(g)]; a = [len(x) for x in S]
    r = [rng.randint(0, a[h]) if rng.random() < .7 else 0 for h in range(g)]
    Rb = [[[rng.randrange(1, P) if lab[v] == h else 0 for v in range(n)] for _ in range(r[h])] for h in range(g)]
    home = [rng.randrange(g) for _ in range(n)]          # alpha(e(omega(v))) : arbitrary
    Phi = []
    for v in range(n):
        vec = [rng.randrange(1, P) if lab[u] == home[v] else 0 for u in range(n)]
        for h in range(g):
            if h != home[v]:
                for bvec in Rb[h]:
                    c = rng.randrange(1, P); vec = [(x + c*y) % P for x, y in zip(vec, bvec)]
        Phi.append(vec)
    X = [home.count(h) for h in range(g)]
    N = sum(max(a[h]-r[h]-X[h], 0) for h in range(g))
    nul = n - rank_mod(Phi); assert nul >= N; tot += 1; eq += nul == N
print(f"generic rows, complete interaction: nullity == N^omega in {eq}/{tot} cases")
