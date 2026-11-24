"""G9c exp2: cost formulas (edge bits U per vertex) for #IS on random d-regular graphs.
Part 1 (asymptotic, per vertex, n -> infinity at fixed d): U/n for
  direct network (= d), Prop. R (Appendix D, all 2^n assignments), Theorem B (class |S| = t <= alpha+1),
  Theorem A (no checks, spectral cut-off), and the enumeration exponent log2 i(G)/n (first moment).
  The same accounting is used for Prop. R and Thm B: U = 2[ nbar*L0 + 2*M*log2(3*sbar) ],
  L0 = max(1, log2(3m/nbar)), beta = h + nbar*L0, M = 3.41*beta, sbar = 16 m/nbar, nbar optimised;
  h = n for Prop. R, h = log2 C(n, alpha) for Thm B.
Part 2 (n = 60 actual graphs): alpha, i(G), lambda, h, certified/heuristic cut-offs.
"""
import math, random, sys
import numpy as np
import networkx as nx
sys.path.insert(0, '/home/user/pnp-recherche/code/gen8/g8b')

H2 = lambda a: 0.0 if a <= 0 or a >= 1 else -(a * math.log2(a) + (1 - a) * math.log2(1 - a))

def fm_exponent(a, d):  # first moment (pairing model), log2 E#IS of density a, per vertex
    if a <= 0: return 0.0
    if a >= 0.5: return -1e9
    return H2(a) + d * ((1 - a) * math.log2(1 - a) - 0.5 * (1 - 2 * a) * math.log2(1 - 2 * a))

def alpha_fm(d):
    lo, hi = 1e-6, 0.4999
    amax = max(np.exp(np.linspace(math.log(1e-9), math.log(0.4999), 40000)), key=lambda a: fm_exponent(a, d))
    lo, hi = amax, 0.4999
    for _ in range(100):
        mid = (lo + hi) / 2
        if fm_exponent(mid, d) > 0: lo = mid
        else: hi = mid
    return lo, fm_exponent(amax, d)

def U_checks(h, mm):  # per vertex: h = class entropy/n, mm = m/n ; returns min over nbar of U/n
    best = None
    for x in np.exp(np.linspace(math.log(1e-6), math.log(mm), 4000)):  # nbar/n
        L0 = max(1.0, math.log2(3 * mm / x))
        beta = h + x * L0
        M = 3.41 * beta
        U = 2 * (x * L0 + 2 * M * math.log2(3 * 16 * mm / x))
        if best is None or U < best[0]: best = (U, x)
    return best

print("Part 1: U/n (edge bits per vertex), random d-regular, asymptotic in n")
print(f"{'d':>6} {'alpha/n':>8} {'enum':>8} {'direct':>7} {'PropR':>8} {'ThmB':>8} {'ThmA':>8} {'B/R':>6} {'enum/ThmA':>10} {'h':>7}")
for d in [4, 8, 16, 32, 64, 128, 256, 1024, 4096, 2**16, 2**20]:
    a, en = alpha_fm(d)
    mm = d / 2
    UR, xr = U_checks(1.0, mm)
    h = H2(a)
    UB, xb = U_checks(h, mm)
    lam = 2 * math.sqrt(d - 1)
    nb = a * (d * a + lam) / 2
    UA = 2 * nb * max(1.0, math.log2(3 * mm / nb))
    print(f"{d:>6} {a:8.5f} {en:8.5f} {d:7d} {UR:8.1f} {UB:8.2f} {UA:8.2f} {UB/UR:6.3f} {en/UA:10.5f} {h:7.4f}")

if len(sys.argv) > 1 and sys.argv[1] == "part1": sys.exit(0)
print("\nPart 2: actual random d-regular graphs, n = 60")
from count_is import to_masks
def is_by_size(n, adj):
    # counts independent sets by size via branching with memo on vertex mask (polynomial vector)
    from functools import lru_cache
    sys.setrecursionlimit(100000)
    @lru_cache(maxsize=None)
    def rec(S):
        if S == 0: return (1,)
        # pick max degree vertex in S
        best, bd = -1, -1
        T = S
        while T:
            b = T & -T; v = b.bit_length() - 1; T ^= b
            dv = (adj[v] & S).bit_count()
            if dv > bd: best, bd = v, dv
        if bd == 0:
            k = S.bit_count()
            return tuple(math.comb(k, j) for j in range(k + 1))
        out = rec(S & ~(1 << best))
        inn = rec(S & ~(1 << best) & ~adj[best])
        res = [0] * max(len(out), len(inn) + 1)
        for j, c in enumerate(out): res[j] += c
        for j, c in enumerate(inn): res[j + 1] += c
        return tuple(res)
    return rec((1 << n) - 1)

def densest_t(G, t, rng, iters=300):
    # heuristic lower bound on max_{|S|=t} e(S): local search swaps
    nodes = list(G.nodes()); best = 0
    for _ in range(6):
        S = set(rng.sample(nodes, t))
        def eS(S): return G.subgraph(S).number_of_edges()
        cur = eS(S)
        for _ in range(iters):
            u = rng.choice(list(S)); v = rng.choice([x for x in nodes if x not in S])
            S2 = (S - {u}) | {v}; e2 = eS(S2)
            if e2 >= cur: S, cur = S2, e2
        best = max(best, cur)
    return best

rng = random.Random(7)
n = 60
print(f"{'d':>3} {'alpha':>5} {'log2 i':>7} {'lam2':>6} {'h=logC(n,a+1)':>14} {'nbar_spec':>9} {'nbar_heur':>9} {'logC(m,nbar_spec)':>17} {'logC(m,nbar_heur)':>17} {'m':>4}")
for d in [4, 6, 8, 12, 16]:
    G = nx.random_regular_graph(d, n, seed=100 + d)
    nn, adj = to_masks(G)
    poly = is_by_size(nn, adj)
    alpha = len(poly) - 1
    while poly[alpha] == 0: alpha -= 1
    tot = sum(poly)
    ev = np.sort(np.linalg.eigvalsh(nx.to_numpy_array(G)))
    lam2 = ev[-2]
    t = alpha + 1
    m = G.number_of_edges()
    nb_spec = int(math.floor(t * (d * t / n + lam2 * (1 - t / n)) / 2))
    nb_heur = densest_t(G, t, rng)
    h = math.log2(math.comb(n, t))
    print(f"{d:>3} {alpha:>5} {math.log2(tot):7.2f} {lam2:6.2f} {h:14.1f} {nb_spec:>9} {nb_heur:>9} {math.log2(math.comb(m, min(nb_spec, m))):17.1f} {math.log2(math.comb(m, nb_heur)):17.1f} {m:>4}")
    sys.stdout.flush()
