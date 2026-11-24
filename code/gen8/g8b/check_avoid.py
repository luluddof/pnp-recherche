# Sanity check of Lemma CA: for any F subset E with |F| < m/(2(1+d t/n)), a t-set minus one endpoint per F-edge
# keeps > nbar edges and avoids F. F chosen as (a) random edges, (b) all edges at a fixed vertex set (stars).
import networkx as nx, random, math
def trial(n, d, nbar, frac, mode, seed):
    rng = random.Random(seed)
    G = nx.random_regular_graph(d, n, seed=seed)
    E = list(G.edges()); m = len(E)
    t = math.ceil(2*math.sqrt(nbar*n/d)) + 3
    bound = m/(2*(1+d*t/n))
    f = int(frac*bound)
    if mode == 'random':
        F = set(map(frozenset, rng.sample(E, f)))
    else:  # stars: edges at vertices 0..; greedy until f edges
        F = set()
        for v in G.nodes():
            for u in G[v]:
                F.add(frozenset((u, v)))
            if len(F) >= f: break
    best = 0
    for _ in range(50):
        X = set(rng.sample(list(G.nodes()), t))
        for e in F:
            u, v = tuple(e)
            if u in X and v in X: X.discard(u)
        ex = sum(1 for u, v in E if u in X and v in X)
        assert not any(frozenset((u, v)) in F for u, v in E if u in X and v in X)
        best = max(best, ex)
    return f, bound, t, best
for (n, d) in [(3000, 8), (3000, 16), (3000, 32)]:
    nbar = int(n/math.log2(2*d))
    for mode in ['random', 'stars']:
        f, b, t, best = trial(n, d, nbar, 0.99, mode, 1)
        print(f"n={n} d={d} nbar={nbar} |F|={f} (bound {b:.0f}, n*sqrt(d log2d)/8={n*math.sqrt(d*math.log2(2*d))/8:.0f}) t={t} mode={mode}: best e(X) avoiding F = {best} > nbar? {best>nbar}")
