import idea, networkx as nx
idea.P = 3
nv=12
G = nx.random_regular_graph(3, nv, seed=1)
eid = {frozenset(e): i + 1 for i, e in enumerate(G.edges())}
ax=[]
for u in G:
    vs = sorted(eid[frozenset((u, w))] for w in G[u])
    ax.append(idea.indicator_wrong_parity(vs, 1 if u in (0,1) else 0))  # even total charge: SAT
for D in [4,5,6]:
    E = idea.pc_closure(ax, len(eid), D)
    print('SAT instance D=%d'%D, 'refuted(BUG)' if frozenset() in E.rows else 'ok not refuted', len(E.rows))
