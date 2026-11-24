"""G9c exp5: Prop. H1 (decision-tree hybrid). L_k(S) = min #leaves of a decision tree on vertex variables of G[S]
all of whose leaves have residual vertex set of size <= k:  L_k(S) = 1 if |S|<=k else min_v L_k(S-v) + L_k(S-N[v]).
Check L_k(V) >= F_Delta(n-k), F(j)=1 (j<=0), F(j)=F(j-1)+F(j-Delta-1); and print growth rates."""
import math, sys, networkx as nx
from functools import lru_cache
sys.setrecursionlimit(1000000)
def F(D, j, memo={}):
    if j <= 0: return 1
    if (D, j) not in memo: memo[(D, j)] = F(D, j - 1) + F(D, j - D - 1)
    return memo[(D, j)]
def minleaves(G, k):
    n = G.number_of_nodes(); adj = [0] * n
    for u, v in G.edges(): adj[u] |= 1 << v; adj[v] |= 1 << u
    @lru_cache(maxsize=None)
    def L(S):
        if S.bit_count() <= k: return 1
        best = None; T = S
        while T:
            b = T & -T; v = b.bit_length() - 1; T ^= b
            c = L(S & ~b) + L(S & ~b & ~adj[v])
            if best is None or c < best: best = c
        return best
    return L((1 << n) - 1)
print(f"{'n':>3} {'d':>2} {'k':>3} {'minleaves':>10} {'F_d(n-k)':>9} ok")
for (n, d) in [(16, 3), (18, 3), (16, 4), (18, 4), (16, 5)]:
    G = nx.random_regular_graph(d, n, seed=n * d)
    for k in [0, n // 4, n // 2]:
        ml = minleaves(G, k); f = F(d, n - k)
        print(f"{n:>3} {d:>2} {k:>3} {ml:>10} {f:>9} {ml >= f}")
        sys.stdout.flush()
print("growth rate log2(x_d)*d/log2(d), x^(d+1)=x^d+1 :")
for d in [4, 16, 64, 256, 1024, 2**16]:
    lo, hi = 1.0, 2.0
    for _ in range(200):
        x = (lo + hi) / 2
        if d * math.log(x) + math.log(x - 1) > 0: hi = x   # x^d (x-1) > 1
        else: lo = x
    print(d, round(math.log2(x) * d / math.log2(d), 4), ' log2 x =', round(math.log2(x), 6))
print("hybrid tree + Prop.R leaves at e*d = 1, units n/d: min_k max{(1-k) W(d)/ln2, k log2(d k/2)}")
for d in [16, 64, 256, 1024, 2**16]:
    lo, hi = 1.0, 2.0
    for _ in range(200):
        x = (lo + hi) / 2
        if d * math.log(x) + math.log(x - 1) > 0: hi = x
        else: lo = x
    best = min(max((1 - k) * d * math.log2(x), k * max(1.0, math.log2(d * k / 2))) for k in [i / 1000 for i in range(1, 1001)])
    print(d, round(best, 3), ' ratio to log2 d:', round(best / math.log2(d), 3))
