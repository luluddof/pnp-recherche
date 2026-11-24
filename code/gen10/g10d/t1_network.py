"""G10d test 1: honest simulation of the Step-3 network (kind-0 marks only) with class profiles s_a,
for independent sets (2-clauses (not x_u or not x_v)).  Value (24) summed over admissible profiles must equal i_t
when nbar >= max_{|S|=t} e(S).  Also checks U bound per profile and that profile-sum == direct identity for
smaller nbar (so the simulation is faithful, including the failure)."""
import itertools, math, random, sys
import networkx as nx
from collections import defaultdict

def brute(G, n):
    edges = list(G.edges())
    it = defaultdict(int); maxe = defaultdict(int); byS = {}
    for mask in range(1 << n):
        t = bin(mask).count('1')
        e = sum(1 for u, v in edges if mask >> u & 1 and mask >> v & 1)
        byS[mask] = e
        if e == 0: it[t] += 1
        maxe[t] = max(maxe[t], e)
    return it, maxe, byS

def identity(byS, n, t, nb):
    tot = 0
    for mask, e in byS.items():
        if bin(mask).count('1') != t: continue
        tot += sum((-1) ** j * math.comb(e, j) for j in range(0, nb + 1))
    return tot

def network_value(G, n, g, part, t, nb):
    """sum over profiles of network values. part[v] in range(g). Returns (Z, maxU over nonzero profiles, bound)."""
    edges = [tuple(sorted(e)) for e in G.edges()]
    m = len(edges)
    # objects = edges; touched nodes sorted; owner = smallest node; test = pair if two distinct nodes
    touch = {e: sorted(set((part[e[0]], part[e[1]]))) for e in edges}
    tests = defaultdict(list)   # pair -> list of edges with that test
    for e in edges:
        if len(touch[e]) == 2: tests[tuple(touch[e])].append(e)
    pairs = [(a, b) for a in range(g) for b in range(a + 1, g)]
    # node tables: key = (s_a, mu_a(owned marked), tuple of mark-strings for each pair involving a) -> sum of W_a
    tables = []
    for a in range(g):
        verts = [v for v in range(n) if part[v] == a]
        eds = [e for e in edges if a in touch[e]]
        T = defaultdict(int)
        for vals in itertools.product([0, 1], repeat=len(verts)):
            x = dict(zip(verts, vals))
            # W_a = 0 if a marked object has a true literal at a: literal (not x_u) true iff x_u = 0, u at a.
            allowed = [e for e in eds if all(x[u] == 1 for u in e if part[u] == a)]
            for r in range(0, min(len(allowed), nb) + 1):
                for marked in itertools.combinations(allowed, r):
                    ms = set(marked)
                    owned = sum(1 for e in ms if touch[e][0] == a)
                    key = []
                    for (p, q) in pairs:
                        if a in (p, q): key.append(tuple(1 if e in ms else 0 for e in tests[(p, q)]))
                        else: key.append(None)
                    T[(sum(vals), owned, tuple(key))] += (-1) ** owned
        tables.append(T)
    # contract: sum over one setting per node with Psi = [marks agree on every test], s sum = t, owned sum <= nb
    # generic g: DP over nodes keeping (s, mu, required strings for pairs with later nodes)
    state = {(0, 0, tuple([None] * len(pairs))): 1}
    profiles_U = 0
    for a in range(g):
        new = defaultdict(int)
        for (s, mu, req), val in state.items():
            for (sa, oa, key), w in tables[a].items():
                if s + sa > t or mu + oa > nb: continue
                ok = True; nreq = list(req)
                for i, (p, q) in enumerate(pairs):
                    if key[i] is None: continue
                    if a == p: nreq[i] = key[i]          # first endpoint fixes the string
                    else:
                        if req[i] != key[i]: ok = False; break
                        nreq[i] = None
                if ok: new[(s + sa, mu + oa, tuple(nreq))] += val * w
        state = new
    Z = sum(v for (s, mu, req), v in state.items() if s == t)
    # U bound: worst admissible profile  U = 2 sum ceil(log C(t_ab, mu_ab)), sum mu_ab <= nb
    # exact max via DP over pairs
    best = {0: 0}
    for pq in pairs:
        tab = len(tests[pq]); nxt = {}
        for used, val in best.items():
            for mu in range(0, min(tab, nb - used) + 1):
                c = math.ceil(math.log2(math.comb(tab, mu))) if math.comb(tab, mu) > 1 else 0
                if nxt.get(used + mu, -1) < val + c: nxt[used + mu] = val + c
        best = nxt
    Umax = 2 * max(best.values())
    bound = 2 * nb * max(1.0, math.log2(3 * m / nb)) + 2 * g * g if nb > 0 else 2 * g * g
    return Z, Umax, bound

rng = random.Random(1)
allok = True
cases = [(10, 3, 3), (12, 3, 3), (12, 4, 3), (12, 3, 4), (14, 3, 3), (14, 4, 4), (12, 5, 3)]
for (n, d, g) in cases:
    for trial in range(2):
        G = nx.random_regular_graph(d, n, seed=rng.randrange(10**6))
        if trial == 1:   # irregular: delete some random edges, add a few
            E = list(G.edges()); rng.shuffle(E); G.remove_edges_from(E[:3])
        perm = list(range(n)); rng.shuffle(perm)
        part = {v: i % g for i, v in enumerate(perm)}
        it, maxe, byS = brute(G, n)
        alpha = max(t for t in it if it[t] > 0)
        for t in range(0, min(alpha + 2, n + 1)):
            nb = maxe[t]
            Z, U, bound = network_value(G, n, g, part, t, nb)
            ok = (Z == it[t]) and (U <= bound)
            # below the cut-off the network must reproduce the (wrong) truncated identity, not i_t
            if nb >= 1:
                Z2, _, _ = network_value(G, n, g, part, t, nb - 1)
                ok2 = (Z2 == identity(byS, n, t, nb - 1)) and (Z2 != it[t])
            else: ok2 = True
            allok &= ok and ok2
            print(f"n={n} d={d} g={g} trial={trial} t={t} i_t={it[t]} nb={nb} Z={Z} U_worst={U} bound={bound:.1f} ok={ok} below_cutoff_wrong={ok2}")
            sys.stdout.flush()
print("ALL OK:", allok)
