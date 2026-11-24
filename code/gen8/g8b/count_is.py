"""Exact #IS counters for benchmarking (G8b).

naive_branch : branch on a max-degree vertex (in: delete N[v]; out: delete v); leaf = edgeless graph.
smart_branch : same + connected components + closed forms for max degree <= 2 + memo on vertex sets.
frontier_dp  : vertex-ordering DP, state = chosen subset of the active frontier (= IS of G[frontier]).
All return (count, work) where work = recursion nodes (branching) or total DP states (frontier).
"""
import random, sys
import networkx as nx
from functools import lru_cache

sys.setrecursionlimit(100000)


def to_masks(G):
    n = G.number_of_nodes()
    idx = {v: i for i, v in enumerate(G.nodes())}
    adj = [0] * n
    for u, v in G.edges():
        adj[idx[u]] |= 1 << idx[v]
        adj[idx[v]] |= 1 << idx[u]
    return n, adj


def bits(m):
    while m:
        b = m & -m
        yield b.bit_length() - 1
        m ^= b


def naive_branch(n, adj):
    work = [0]

    def rec(S):
        work[0] += 1
        best, bd = -1, 0
        for v in bits(S):
            dv = bin(adj[v] & S).count("1")
            if dv > bd:
                best, bd = v, dv
        if bd == 0:
            return 1 << bin(S).count("1")
        out = rec(S & ~(1 << best))
        inn = rec(S & ~(1 << best) & ~adj[best])
        return out + inn

    c = rec((1 << n) - 1)
    return c, work[0]


# closed forms: #IS of path on k vertices = F(k+2); cycle on k vertices = L(k)
def fib(k):
    a, b = 0, 1
    for _ in range(k):
        a, b = b, a + b
    return a


def lucas(k):
    a, b = 2, 1
    for _ in range(k):
        a, b = b, a + b
    return a


def smart_branch(n, adj):
    work = [0]
    memo = {}

    def components(S):
        comps = []
        rem = S
        while rem:
            v = (rem & -rem).bit_length() - 1
            comp = 1 << v
            frontier = comp
            while frontier:
                nb = 0
                for u in bits(frontier):
                    nb |= adj[u]
                nb &= rem & ~comp
                comp |= nb
                frontier = nb
            comps.append(comp)
            rem &= ~comp
        return comps

    def rec(S):
        if S in memo:
            return memo[S]
        work[0] += 1
        res = 1
        for C in components(S):
            res *= solve_conn(C)
        memo[S] = res
        return res

    def solve_conn(C):
        if C in memo:
            return memo[C]
        k = bin(C).count("1")
        best, bd = -1, -1
        for v in bits(C):
            dv = bin(adj[v] & C).count("1")
            if dv > bd:
                best, bd = v, dv
        if bd <= 2:
            m = sum(bin(adj[v] & C).count("1") for v in bits(C)) // 2
            r = lucas(k) if m == k and k >= 3 else fib(k + 2)
            memo[C] = r
            return r
        work[0] += 1
        r = rec(C & ~(1 << best)) + rec(C & ~(1 << best) & ~adj[best])
        memo[C] = r
        return r

    c = rec((1 << n) - 1)
    return c, work[0]


def greedy_order(n, adj):
    """order minimizing the number of active vertices greedily (pathwidth heuristic)."""
    processed = 0
    order = []
    for _ in range(n):
        best, bkey = None, None
        for v in range(n):
            if processed >> v & 1:
                continue
            P = processed | (1 << v)
            # active = processed vertices with an unprocessed neighbour
            act = 0
            for u in bits(P):
                if adj[u] & ~P:
                    act += 1
            key = (act, -bin(adj[v] & processed).count("1"))
            if bkey is None or key < bkey:
                best, bkey = v, key
        order.append(best)
        processed |= 1 << best
    return order


def frontier_dp(n, adj, order=None):
    if order is None:
        order = greedy_order(n, adj)
    states = {0: 1}  # chosen subset of active vertices -> count
    processed = 0
    work = 0
    maxst = 1
    for v in order:
        processed |= 1 << v
        new = {}
        for s, c in states.items():
            new[s] = new.get(s, 0) + c
            if not (adj[v] & s):
                t = s | (1 << v)
                new[t] = new.get(t, 0) + c
        # retire vertices with no unprocessed neighbours
        retire = 0
        for u in bits(processed):
            if not (adj[u] & ~processed):
                retire |= 1 << u
        if retire:
            red = {}
            for s, c in new.items():
                t = s & ~retire
                red[t] = red.get(t, 0) + c
            new = red
        processed_keep = processed
        states = new
        work += len(states)
        maxst = max(maxst, len(states))
    return sum(states.values()), work, maxst


def brute(n, adj):
    c = 0
    for S in range(1 << n):
        ok = True
        for v in bits(S):
            if adj[v] & S:
                ok = False
                break
        c += ok
    return c


if __name__ == "__main__":
    for seed in range(20):
        rng = random.Random(seed)
        n = rng.randint(4, 16)
        d = rng.choice([2, 3, 4, 5])
        if n * d % 2:
            n += 1
        if d >= n:
            continue
        G = nx.random_regular_graph(d, n, seed=seed)
        nn, adj = to_masks(G)
        b = brute(nn, adj)
        a1 = naive_branch(nn, adj)[0]
        a2 = smart_branch(nn, adj)[0]
        a3 = frontier_dp(nn, adj)[0]
        assert b == a1 == a2 == a3, (n, d, b, a1, a2, a3)
    print("all counters agree with brute force")
