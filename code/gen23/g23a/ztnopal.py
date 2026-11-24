import sys
from zt import *; from g64 import su3_sylow, gf
from pysat.solvers import Cadical195
q = int(sys.argv[1]); G = su3_sylow(q); mul, inv, el = G; N = len(mul)
n = {2: 2, 4: 4, 8: 6}[q]; poly = {2: 0b111, 4: 0b10011, 8: 0b1000011}[q]; m = gf(n, poly)
Z = frozenset(x for x in range(N) if el[x][0] == 0); Zs = sorted(Z)
Z0 = next(c for c in central_subgroups(mul, inv, len(Z) // 2) if c <= Z)
cos = {}; reps = []
for x in range(N):
    k = min(mul[x][z] for z in Zs)
    if k not in cos: cos[k] = len(reps); reps.append(k)
for mu in range(2, 1 << n):
    allowed = [[x for x in range(N) if el[x][0] == m(mu, el[reps[c]][0])] for c in range(len(reps))]
    E = ZT(mul, inv, Z, Z0, allowed, nopal=True); S = Cadical195(bootstrap_with=E.cl)
    r = S.solve(); print('mu', mu, 'colour structure without palette condition:', r, flush=True)
    if r:
        tb, lam = E.decode(S.get_model()); col = E.array(tb, lam)
        # intercalation graph bipartite?  colours = (block, idx); check which block pairs are adjacent
        import itertools
        adj = set()
        cls = {}
        for a in range(N):
            for b in range(N): cls.setdefault(col[a][b], []).append((a, b))
        for k, lst in cls.items():
            for (a, b), (a2, b2) in itertools.combinations(lst, 2):
                k2 = col[a][b2]; adj.add((k[1] // N, k2[1] // N))
        print('   block adjacency (block of colour, block of partner):', sorted(adj), 'sign ok:', sign_solve([[(x[1],) for x in row] for row in col]) is not None)
        break
