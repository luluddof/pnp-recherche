# Exact treewidth of the line graph L(K_K) (primal graph of the edge variables of a (K,lam)-clique) by the
# O*(2^n) DP of Bodlaender et al.:  TW(S) = min_{v in S} max(TW(S-v), |Q(S-v, v)|).
# Variable elimination with any order enumerates >= tw+1 edge strings at once; compare with s(K) = floor((K+1)^2/4)-1.
import itertools, sys
from functools import lru_cache
def line_graph(K):
    E = list(itertools.combinations(range(K), 2)); n = len(E)
    adj = [0] * n
    for i in range(n):
        for j in range(n):
            if i != j and set(E[i]) & set(E[j]): adj[i] |= 1 << j
    return n, adj
def Q(S, v, adj, n):
    # vertices outside S and != v reachable from v via paths whose interior lies in S
    seen = 1 << v; stack = [v]; out = 0
    while stack:
        u = stack.pop(); nb = adj[u] & ~seen
        seen |= nb
        inS = nb & S; out |= nb & ~S
        x = inS
        while x:
            b = x & -x; stack.append(b.bit_length() - 1); x ^= b
    return bin(out & ~(1 << v)).count("1")
def tw(K):
    n, adj = line_graph(K)
    TW = {0: -1}
    for size in range(1, n + 1):
        for comb in itertools.combinations(range(n), size):
            S = 0
            for c in comb: S |= 1 << c
            best = 10**9
            for v in comb:
                Sm = S & ~(1 << v)
                val = max(TW[Sm], Q(Sm, v, adj, n))
                if val < best: best = val
            TW[S] = best
    return TW[(1 << n) - 1]
for K in range(3, 7):
    t = tw(K); print(f"K={K}: tw(L(K_K))={t}  tw+1={t+1}  s(K)={(K+1)**2//4-1}", flush=True)
