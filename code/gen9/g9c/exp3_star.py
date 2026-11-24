"""G9c exp3: star identity  prod_{v in N+(u)} (1 - x_u x_v) = 1 - x_u + x_u prod_{v in N+(u)} (1 - x_v).
Edge bits of the resulting network = (vertex, foreign cluster) incidences:
  inc  = sum_u #{clusters b != c(u) : N(u)  meets A_b}   (undirected stars)
  inc+ = sum_u #{clusters b != c(u) : N+(u) meets A_b}   (each edge owned by its smaller endpoint)
Compare: random balanced clustering, a local-search clustering, the formula n (K-1)(1-(1-1/K)^d),
and the spectral lower bound  n d^2/(d^2 smax/n + lam^2) - n  (Cauchy-Schwarz + trace of A^2 on clusters).
"""
import math, random, sys
import numpy as np, networkx as nx

def incidences(adj, lab, K, oriented=False):
    tot = 0
    for u, nb in enumerate(adj):
        s = set(lab[v] for v in nb if (not oriented or v > u))
        s.discard(lab[u])
        tot += len(s)
    return tot

def local_search(adj, lab, K, rng, sweeps=6):
    n = len(adj); lab = lab[:]
    def inc_local(vs):
        t = 0
        for u in vs:
            s = set(lab[v] for v in adj[u]); s.discard(lab[u]); t += len(s)
        return t
    for _ in range(sweeps * n):
        u, w = rng.randrange(n), rng.randrange(n)
        if lab[u] == lab[w]: continue
        aff = set([u, w]) | set(adj[u]) | set(adj[w])
        before = inc_local(aff)
        lab[u], lab[w] = lab[w], lab[u]
        if inc_local(aff) > before:
            lab[u], lab[w] = lab[w], lab[u]
    return lab

rng = random.Random(5)
n = 2048
print(f"n={n}; all numbers divided by n*d (so 1.0 = one bit per edge-endpoint, direct network = 1.0)")
print(f"{'d':>3} {'K':>5} {'rand inc':>9} {'formula':>8} {'search inc':>10} {'rand inc+':>9} {'spectral LB':>11} {'lam':>6}")
for d in [8, 16, 32]:
    G = nx.random_regular_graph(d, n, seed=d)
    adj = [list(G[u]) for u in range(n)]
    ev = np.linalg.eigvalsh(nx.to_numpy_array(G)); lam = max(abs(ev[0]), abs(ev[-2]))
    for K in [2, 4, 8, 16, 32, 64, 256]:
        lab = [i % K for i in range(n)]; rng.shuffle(lab)
        inc = incidences(adj, lab, K); incp = incidences(adj, lab, K, True)
        form = n * (K - 1) * (1 - (1 - 1 / K) ** d)
        lab2 = local_search(adj, lab, K, rng) if K <= 64 else lab
        inc2 = incidences(adj, lab2, K)
        smax = n / K
        lb = max(0.0, n * d * d / (d * d * smax / n + lam * lam) - n)
        print(f"{d:>3} {K:>5} {inc/(n*d):9.3f} {form/(n*d):8.3f} {inc2/(n*d):10.3f} {incp/(n*d):9.3f} {lb/(n*d):11.3f} {lam:6.2f}")
        sys.stdout.flush()
