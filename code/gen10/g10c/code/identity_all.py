# G10c independent check of the counting identity (definitions only).
# For a graph G, terminals K, l, modulus p:
#   S_ST  = sum over (Y,X): K<=Y, X<=E(G[Y]), |Y|=l+1, |X|=l  of p^{cc(Y,X)}
#   claim: S_ST = p*|ST(l)| mod p^2, where |ST(l)| is computed INDEPENDENTLY by Kirchhoff:
#          |ST(l)| = sum_{Y>=K, |Y|=l+1} tau(G[Y]).
# cycles: (Y,X) 2-regular on Y, |Y|=l ; independent count by DFS.
# paths : |Y|=l, |X|=l-1, all degrees<=2 ; independent count by DFS.
import itertools, sys
from fractions import Fraction
import networkx as nx
from networkx.generators.atlas import graph_atlas_g

def cc_of(Ymask, Xedges, n):
    par = list(range(n))
    def f(a):
        while par[a] != a:
            par[a] = par[par[a]]; a = par[a]
        return a
    for u, v in Xedges:
        par[f(u)] = f(v)
    return len({f(v) for v in range(n) if Ymask >> v & 1})

def det_int(M):
    k = len(M)
    if k == 0: return 1
    A = [[Fraction(x) for x in r] for r in M]; d = Fraction(1)
    for i in range(k):
        piv = next((r for r in range(i, k) if A[r][i] != 0), None)
        if piv is None: return 0
        if piv != i: A[i], A[piv] = A[piv], A[i]; d = -d
        d *= A[i][i]
        for r in range(i + 1, k):
            fct = A[r][i] / A[i][i]
            if fct:
                for c in range(i, k): A[r][c] -= fct * A[i][c]
    assert d.denominator == 1
    return int(d)

def tau(Y, E, wt=None):
    Y = list(Y); idx = {v: i for i, v in enumerate(Y)}; k = len(Y)
    L = [[0] * k for _ in range(k)]
    for j, (u, v) in enumerate(E):
        if u in idx and v in idx:
            w = 1 if wt is None else wt[j]
            a, b = idx[u], idx[v]
            L[a][a] += w; L[b][b] += w; L[a][b] -= w; L[b][a] -= w
    return det_int([r[1:] for r in L[1:]])

def count_cycles(n, adj, l):
    # undirected simple cycles with exactly l vertices: start at min vertex, each counted twice
    cnt = 0
    def dfs(start, v, used, depth):
        nonlocal cnt
        if depth == l:
            if start in adj[v]: cnt += 1
            return
        for w in adj[v]:
            if w > start and not used >> w & 1:
                dfs(start, w, used | 1 << w, depth + 1)
    for s in range(n): dfs(s, s, 1 << s, 1)
    return cnt // 2 if l >= 3 else 0

def count_paths(n, adj, l):
    # undirected simple paths on exactly l vertices (l>=2), each counted twice by DFS
    cnt = 0
    def dfs(v, used, depth):
        nonlocal cnt
        if depth == l: cnt += 1; return
        for w in adj[v]:
            if not used >> w & 1: dfs(w, used | 1 << w, depth + 1)
    for s in range(n): dfs(s, 1 << s, 1)
    return cnt // 2

def check_graph(n, E, primes=(2, 3, 5, 7), allK=True, Ks=None):
    m = len(E)
    adj = [set() for _ in range(n)]
    for u, v in E: adj[u].add(v); adj[v].add(u)
    # tables: A[Y][l] = dict cc->count over X<=E(G[Y]) with |X|=l (no degree constraint)
    #         Cyc[l] = dict cc->count over 2-regular (Y,X) with |Y|=l
    #         Pth[l] = dict cc->count over (Y,X), |Y|=l,|X|=l-1, maxdeg<=2
    A = {}; Cyc = {}; Pth = {}
    for Y in range(1, 1 << n):
        EY = [e for e in E if Y >> e[0] & 1 and Y >> e[1] & 1]
        ny = bin(Y).count('1')
        for r in range(0, len(EY) + 1):
            if r not in (ny - 1, ny): continue   # only sizes that are used
            for X in itertools.combinations(EY, r):
                c = cc_of(Y, X, n)
                deg = [0] * n
                for u, v in X: deg[u] += 1; deg[v] += 1
                if r == ny - 1:
                    d = A.setdefault(Y, {}); d[c] = d.get(c, 0) + 1
                    if max(deg) <= 2 or r == 0:
                        d = Pth.setdefault(ny, {}); d[c] = d.get(c, 0) + 1
                if r == ny and all(deg[v] == 2 for v in range(n) if Y >> v & 1):
                    d = Cyc.setdefault(ny, {}); d[c] = d.get(c, 0) + 1
    tests = 0
    if Ks is None: Ks = range(1, 1 << n) if allK else [1]
    tauY = {Y: tau([v for v in range(n) if Y >> v & 1], E) for Y in range(1, 1 << n)}
    for K in Ks:
        for l in range(1, n):
            Ys = [Y for Y in range(1, 1 << n) if Y & K == K and bin(Y).count('1') == l + 1]
            trees = sum(tauY[Y] for Y in Ys)
            mincc = min([c for Y in Ys for c in A.get(Y, {})], default=None)
            assert mincc is None or mincc >= 1
            for p in primes + (4, 6):          # also composite moduli: Theorem 1 says "every integer p>=2"
                S = sum(cnt * p ** c for Y in Ys for c, cnt in A.get(Y, {}).items())
                assert S % p == 0 and (S // p) % p == trees % p, ("ST", n, E, K, l, p, S, trees)
                tests += 1
    for l in range(3, n + 1):
        cyc = count_cycles(n, adj, l)
        for p in primes + (4, 6):
            S = sum(cnt * p ** c for c, cnt in Cyc.get(l, {}).items())
            assert all(c >= 1 for c in Cyc.get(l, {})) and Cyc.get(l, {}).get(1, 0) == cyc
            assert S % p == 0 and (S // p) % p == cyc % p, ("cyc", n, E, l, p); tests += 1
    for l in range(2, n + 1):
        pth = count_paths(n, adj, l)
        for p in primes + (4, 6):
            S = sum(cnt * p ** c for c, cnt in Pth.get(l, {}).items())
            assert Pth.get(l, {}).get(1, 0) == pth, ("pathcount", n, E, l, Pth.get(l), pth)
            assert S % p == 0 and (S // p) % p == pth % p, ("path", n, E, l, p); tests += 1
    return tests

if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "atlas":
        tot = 0; ng = 0
        for G in graph_atlas_g():
            n = G.number_of_nodes()
            if n < 2 or n > 6: continue
            tot += check_graph(n, sorted(tuple(sorted(e)) for e in G.edges())); ng += 1
        print(f"atlas: all {ng} non-isomorphic graphs with 2..6 vertices, all terminal sets K, all l: {tot} identity checks passed")
    else:
        import random
        rnd = random.Random(20261124); tot = 0; ng = 0
        for n in (7, 8, 9):
            for trial in range(int(sys.argv[2])):
                m = rnd.randint(n - 1, min(n * (n - 1) // 2, 14))
                allE = list(itertools.combinations(range(n), 2)); rnd.shuffle(allE); E = sorted(allE[:m])
                # random terminal sets: relabel so that K varies -> use 6 random K
                Ks = None if n <= 7 else [rnd.randrange(1, 1 << n) for _ in range(12)]
                tot += check_graph(n, E, Ks=Ks); ng += 1
        print(f"random: {ng} random graphs n=7..9, m<=14: {tot} identity checks passed")
