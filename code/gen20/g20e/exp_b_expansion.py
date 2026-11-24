"""Exp B: edge expansion of the graphs the two algorithms contract, for random 3-CNF at density D.
Graph G_chain (Claude §2.1): nodes = clauses, for each variable its occurrences joined in a chain (random order) -> M = sum(d_v - 1) wires.
Graph G_primal: variables, edge per co-occurrence (multigraph weights).
Reported, as a fraction of the number of wires M:
  spec  = lambda_2(L) * n_nodes / 4 / M      (lower bound on the min BISECTION width / M; [Fiedler/Donath-Hoffman])
  kl    = best of several Kernighan-Lin bisections / M (upper bound on min bisection / M)
  cheeger lower bound on conductance phi >= lambda_2(normalised L)/2.
Also: after deleting the b highest-degree variables (branching on them; their wires vanish), the same numbers.
"""
import numpy as np, random, sys, networkx as nx, scipy.sparse as sp, scipy.sparse.linalg as sla

def rand_cnf(n, m, rng):
    return [tuple((v, rng.random() < 0.5) for v in rng.sample(range(n), 3)) for _ in range(m)]

def chain_graph(n, cls, removed=frozenset(), rng=None):
    occ = {v: [] for v in range(n)}
    for j, c in enumerate(cls):
        for v, _ in c:
            if v not in removed: occ[v].append(j)
    G = nx.MultiGraph(); G.add_nodes_from(range(len(cls)))
    for v, l in occ.items():
        rng.shuffle(l)
        for a, b in zip(l, l[1:]): G.add_edge(a, b)
    return G

def simple_weighted(G):
    H = nx.Graph(); H.add_nodes_from(G.nodes())
    for a, b in G.edges():
        if H.has_edge(a, b): H[a][b]['weight'] += 1
        else: H.add_edge(a, b, weight=1)
    return H

def measures(G, rng, tries=6):
    H = simple_weighted(G)
    H.remove_nodes_from([v for v in list(H.nodes()) if H.degree(v) == 0])
    M = H.size(weight='weight'); N = H.number_of_nodes()
    if not nx.is_connected(H):
        comp = max(nx.connected_components(H), key=len)
        note = f'disconnected: giant {len(comp)}/{N}'
    else: note = ''
    L = nx.laplacian_matrix(H, weight='weight').astype(float)
    vals = sla.eigsh(L, k=2, sigma=-1e-3, which='LM', return_eigenvectors=False)
    l2 = sorted(vals)[1]
    NL = nx.normalized_laplacian_matrix(H, weight='weight').astype(float)
    nvals = sla.eigsh(NL, k=2, sigma=-1e-3, which='LM', return_eigenvectors=False)
    nl2 = sorted(nvals)[1]
    best = None
    for t in range(tries):
        a, b = nx.algorithms.community.kernighan_lin_bisection(H, weight='weight', seed=rng.randrange(10 ** 9))
        cut = nx.cut_size(H, a, b, weight='weight')
        best = cut if best is None else min(best, cut)
    return dict(M=M, N=N, spec=l2 * N / 4 / M, kl=best / M, cheeger=nl2 / 2, note=note)

if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
    D = 4.26
    print('graph n M  spec_lb(bisect/M)  KL_ub(bisect/M)  cheeger_lb(phi)  [b removed]')
    for n in [100, 200, 400]:
        cls = rand_cnf(n, int(round(D * n)), rng)
        deg = np.zeros(n, int)
        for c in cls:
            for v, _ in c: deg[v] += 1
        order = list(np.argsort(-deg))
        for frac in [0.0, 0.05, 0.1, 0.2, 0.3]:
            b = int(frac * n)
            rem = frozenset(order[:b])
            G = chain_graph(n, cls, rem, rng)
            r = measures(G, rng, tries=4 if n >= 400 else 6)
            print('chain', n, int(r['M']), f"{r['spec']:.4f}", f"{r['kl']:.4f}", f"{r['cheeger']:.4f}", f'b={b}', r['note'], flush=True)
        P = nx.MultiGraph(); P.add_nodes_from(range(n))
        for c in cls:
            vs = [v for v, _ in c]
            P.add_edges_from([(vs[0], vs[1]), (vs[0], vs[2]), (vs[1], vs[2])])
        r = measures(P, rng)
        print('primal', n, int(r['M']), f"{r['spec']:.4f}", f"{r['kl']:.4f}", f"{r['cheeger']:.4f}", r['note'], flush=True)
