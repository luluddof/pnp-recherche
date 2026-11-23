"""G6a: clause-product interpolation schemes (Limits paper, Def. 5.1) -- exact checker.

A design is: formula (clauses), points X (nonzero ints), degree R, for each occurrence o=(e,v)
a polynomial vector P_o(T) (dict u -> list of Fraction coeffs, degree <= R, P_o(0)=e_v),
for each t: partition groups[t] (list of lists), alpha[t][e] (group index of clause e),
and measurement spaces Rsp[t][h] (list of row vectors over S_{t,h}, Fractions).
"""
import random, itertools, math
from fractions import Fraction as Fr

# ---------------------------------------------------------------- formulas
def random_3regular_formula(n, seed):
    """n variables, n clauses, every variable in exactly 3 clauses, clauses have 3 distinct vars."""
    rng = random.Random(seed)
    while True:
        stubs = [v for v in range(n) for _ in range(3)]
        rng.shuffle(stubs)
        cl = [stubs[3*i:3*i+3] for i in range(n)]
        if all(len(set(c)) == 3 for c in cl):
            return [[(v, rng.random() < 0.5) for v in c] for c in cl]  # (var, positive?)

def occurrences(F):
    return [(e, v) for e, c in enumerate(F) for (v, _) in c]

# ---------------------------------------------------------------- polynomials
def pev(c, t):
    r = Fr(0)
    for a in reversed(c):
        r = r * t + a
    return r

def Pval(P, t):
    """P: dict u -> coeff list. returns dict u -> value (nonzero only)."""
    out = {}
    for u, c in P.items():
        x = pev(c, t)
        if x != 0:
            out[u] = x
    return out

# ---------------------------------------------------------------- linear algebra over Q
def rank_Q(rows):
    M = [list(r) for r in rows if any(x != 0 for x in r)]
    if not M:
        return 0
    m = len(M[0]); rk = 0; col = 0
    while rk < len(M) and col < m:
        piv = next((i for i in range(rk, len(M)) if M[i][col] != 0), None)
        if piv is None:
            col += 1; continue
        M[rk], M[piv] = M[piv], M[rk]
        p = M[rk][col]
        for i in range(rk + 1, len(M)):
            if M[i][col] != 0:
                f = M[i][col] / p
                M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
        rk += 1; col += 1
    return rk

def in_span(basis, vec):
    return rank_Q(basis + [vec]) == rank_Q(basis)

# ---------------------------------------------------------------- design checks and costs
def check_design(F, n, X, R, P, groups, alpha, Rsp):
    """Verify P_o(0)=e_v, deg<=R, and condition (3) exactly. Return per-point stats."""
    occ = occurrences(F)
    for o in occ:
        e, v = o
        for u, c in P[o].items():
            assert len(c) <= R + 1
            assert c[0] == (1 if u == v else 0), (o, u, c)
        assert v in P[o] and P[o][v][0] == 1
    stats = []
    for t in X:
        gid = {}
        for h, S in enumerate(groups[t]):
            for u in S:
                gid[u] = h
        assert len(gid) == n
        for o in occ:
            e, v = o
            val = Pval(P[o], t)
            for h, S in enumerate(groups[t]):
                if h == alpha[t][e]:
                    continue
                vec = [val.get(u, Fr(0)) for u in S]
                if any(x != 0 for x in vec):
                    assert in_span(Rsp[t][h], vec), ("condition (3) fails", t, o, h)
        rho = sum(rank_Q(Rsp[t][h]) for h in range(len(groups[t])))
        s = max(len(S) for S in groups[t])
        mu = sum(1 for h in range(len(groups[t])) if rank_Q(Rsp[t][h]) > 0)
        stats.append(dict(t=t, s=s, rho=rho, mu=mu, g=len(groups[t])))
    return stats

def interaction_width(F, X, P, groups, alpha):
    import networkx as nx
    from networkx.algorithms.approximation import treewidth_min_fill_in
    ws = []
    for t in X:
        G = nx.Graph(); G.add_nodes_from(range(len(groups[t])))
        for e, c in enumerate(F):
            cl = {alpha[t][e]}
            for (v, _) in c:
                val = Pval(P[(e, v)], t)
                for h, S in enumerate(groups[t]):
                    if h != alpha[t][e] and any(u in val for u in S):
                        cl.add(h)
            cl = list(cl)
            G.add_edges_from(itertools.combinations(cl, 2))
        ws.append(treewidth_min_fill_in(G)[0])
    return ws

# ---------------------------------------------------------------- Theorem B quantities
def expected_deficiency(F, n, groups_t, alpha_t):
    """E_omega sum_h (|S_h| - |I^omega_h|)_+ exactly (float), plus sigma' and mean-imbalance term."""
    occ_of = {v: [] for v in range(n)}
    for e, c in enumerate(F):
        for (v, _) in c:
            occ_of[v].append(e)
    g = len(groups_t)
    pis = {}
    split = 0
    for v in range(n):
        d = len(occ_of[v])
        cnt = {}
        for e in occ_of[v]:
            cnt[alpha_t[e]] = cnt.get(alpha_t[e], 0) + 1
        if len(cnt) >= 2:
            split += 1
        pis[v] = {h: k / d for h, k in cnt.items()}
    Edef = 0.0; imb = 0.0
    for h in range(g):
        a = len(groups_t[h])
        dist = [1.0]
        for v in range(n):
            p = pis[v].get(h, 0.0)
            if p == 0.0:
                continue
            nd = [0.0] * (len(dist) + 1)
            for k, q in enumerate(dist):
                nd[k] += q * (1 - p); nd[k + 1] += q * p
            dist = nd
        Edef += sum(q * max(a - k, 0) for k, q in enumerate(dist))
        mu = sum(k * q for k, q in enumerate(dist))
        imb += abs(mu - a)
    return Edef, split, imb / 2

def phi_matrix_nullity(F, n, t, P, omega):
    rows = []
    for v in range(n):
        val = Pval(P[omega[v]], t)
        rows.append([val.get(u, Fr(0)) for u in range(n)])
    return n - rank_Q(rows)

def deficiency(groups_t, alpha_t, omega, n):
    cnt = [0] * len(groups_t)
    for v in range(n):
        cnt[alpha_t[omega[v][0]]] += 1
    return sum(max(len(S) - cnt[h], 0) for h, S in enumerate(groups_t))
