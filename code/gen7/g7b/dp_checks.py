"""G7b: DP-type cost model checks. Imports G6a's exact checker (read-only)."""
import sys, math, random, itertools
sys.path.insert(0, "/home/user/pnp-recherche/code/gen6/g6a")
from fractions import Fraction as Fr
from scheme import (random_3regular_formula, occurrences, Pval, rank_Q, check_design,
                    phi_matrix_nullity)
from designs import (design_cut_measure, design_global_sum, design_random_dense)

C0 = 0.216

# ----------------------------------------------------------- exact (weighted) treewidth
def weighted_tw(nodes, adj, w):
    """min over elimination orders of max_v (w(v) + w(Q(S,v))) = min over tree decompositions
    of the max bag weight (Bodlaender et al. subset DP, exact). Returns max bag weight."""
    nodes = list(nodes); k = len(nodes)
    if k == 0:
        return 0
    idx = {v: i for i, v in enumerate(nodes)}
    nb = [0] * k
    for v in nodes:
        for u in adj[v]:
            if u in idx and u != v:
                nb[idx[v]] |= 1 << idx[u]
    W = [w[v] for v in nodes]
    def Q(S, i):
        # vertices outside S∪{i} reachable from i via paths whose interior lies in S
        seen = 1 << i; stack = [i]; out = 0
        while stack:
            j = stack.pop()
            m = nb[j] & ~seen
            while m:
                b = m & -m; m ^= b; seen |= b
                jj = b.bit_length() - 1
                if S >> jj & 1:
                    stack.append(jj)
                else:
                    out |= b
        return out
    INF = float("inf")
    best = [INF] * (1 << k); best[0] = 0
    for S in range(1, 1 << k):
        m = S
        while m:
            b = m & -m; m ^= b; i = b.bit_length() - 1
            R = S ^ b
            if best[R] >= best[S]:
                continue
            q = Q(R, i)
            cost = W[i] + sum(W[j] for j in range(k) if q >> j & 1)
            c = max(best[R], cost)
            if c < best[S]:
                best[S] = c
    return best[(1 << k) - 1]

def primal_adj(F, n):
    adj = {v: set() for v in range(n)}
    for c in F:
        vs = [v for v, _ in c]
        for a, b in itertools.permutations(vs, 2):
            adj[a].add(b)
    return adj

def mmd_plus_lb(adj):
    """minor-min-width lower bound on treewidth (contract min-degree vertex into min-degree nbr)."""
    G = {v: set(s) for v, s in adj.items()}
    lb = 0
    while len(G) > 1:
        v = min(G, key=lambda x: len(G[x]))
        lb = max(lb, len(G[v]))
        if not G[v]:
            del G[v]; continue
        u = min(G[v], key=lambda x: len(G[x]))
        for x in G[v]:
            if x != u:
                G[x].discard(v); G[x].add(u); G[u].add(x)
        G[u].discard(v); G[u].discard(u); del G[v]
    return lb

def interaction_graph(F, t, P, groups_t, alpha_t):
    gid = {u: h for h, S in enumerate(groups_t) for u in S}
    adj = {h: set() for h in range(len(groups_t))}
    for e, c in enumerate(F):
        cl = {alpha_t[e]}
        for (v, _) in c:
            for u in Pval(P[(e, v)], t):
                cl.add(gid[u])
        for a, b in itertools.permutations(cl, 2):
            adj[a].add(b)
    return adj

def torso(adj, measured):
    T = {h: set(x for x in adj[h] if x in measured) for h in measured}
    for h in adj:
        if h not in measured:
            N = [x for x in adj[h] if x in measured]
            assert all(x in measured for x in adj[h]), "two unmeasured groups adjacent"
            for a, b in itertools.permutations(N, 2):
                T[a].add(b)
    return T

# ----------------------------------------------------------- Part 1: literal DP model, trivial design
def trivial_design(F, n, R, X):
    s = math.isqrt(n - 1) + 1          # ceil(sqrt n)
    g = -(-n // s)
    groups_t = [list(range(h * s, min(n, (h + 1) * s))) for h in range(g)]
    gid = {u: h for h, S in enumerate(groups_t) for u in S}
    alpha_t = [gid[c[0][0]] for c in F]
    P = {o: {o[1]: [Fr(1)] + [Fr(0)] * R} for o in occurrences(F)}   # constant curve e_v
    groups = {t: groups_t for t in X}; alpha = {t: alpha_t for t in X}
    Rsp = {t: [[[Fr(int(i == j)) for j in range(len(S))] for i in range(len(S))] for S in groups_t] for t in X}
    return P, groups, alpha, Rsp

def part1():
    print("== Part 1: literal DP-type model. Trivial design (P_o = e_v, full measurement, q_t = 0)")
    for n in [12, 16, 20, 24, 30]:
        F = random_3regular_formula(n, seed=100 + n)
        R = 1; M = 3 * n
        X = list(range(1, M * R + 2))[:3] if False else [1]   # checks at one representative point
        P, groups, alpha, Rsp = trivial_design(F, n, R, list(range(1, M * R + 2)))
        Xall = list(range(1, M * R + 2))                      # |X| = MR+1, q_t = 0: Lagrange exact
        st = check_design(F, n, [1, 2, Xall[-1]], R, P, groups, alpha, Rsp)   # (3) exact (design is t-independent)
        t = 1
        I = interaction_graph(F, t, P, groups[t], alpha[t])
        w = weighted_tw(I.keys(), I, {h: 1 for h in I}) - 1
        meas = set(I)
        Tor = torso(I, meas)
        W = weighted_tw(Tor.keys(), Tor, {h: len(groups[t][h]) for h in Tor})  # reading-weighted width
        adj = primal_adj(F, n)
        lb = mmd_plus_lb(adj)
        # exact tw(Pr F) only for small n
        tau = weighted_tw(range(n), adj, {v: 1 for v in range(n)}) - 1 if n <= 16 else None
        th81 = math.sqrt((lb + 1) / 2) - 1
        Zt = sum(1 for o in occurrences(F) if Pval(P[o], t).get(o[1], 0) == 0)
        print(f"  n={n:2d} |X|={len(Xall)} q_t=0 |Z_t|={Zt}  s={st[0]['s']} g={st[0]['g']} rho={st[0]['rho']}"
              f"  w_t=tw(I_t)={w}  literal exponent max(s,w)={max(st[0]['s'], w)}  ceil(sqrt n)={math.isqrt(n-1)+1}"
              f"  | tw(Pr F): exact={tau} mmd+lb={lb}; Thm 8.1 exponent >= {th81:.2f}"
              f"  | reading-weighted W_t={W} (>= tw(PrF)+1-s: {W >= (tau if tau is not None else lb) + 1 - st[0]['s']})")

# ----------------------------------------------------------- Part 2: single-point lemma, per-group bound
def group_dists(F, n, groups_t, alpha_t):
    occ_of = {v: [] for v in range(n)}
    for e, c in enumerate(F):
        for (v, _) in c:
            occ_of[v].append(e)
    pis = {}
    for v in range(n):
        d = len(occ_of[v]); cnt = {}
        for e in occ_of[v]:
            cnt[alpha_t[e]] = cnt.get(alpha_t[e], 0) + 1
        pis[v] = {h: k / d for h, k in cnt.items()}
    out = []
    for h in range(len(groups_t)):
        dist = [1.0]; V = 0.0
        for v in range(n):
            p = pis[v].get(h, 0.0)
            if p == 0.0:
                continue
            V += p * (1 - p)
            nd = [0.0] * (len(dist) + 1)
            for k, q in enumerate(dist):
                nd[k] += q * (1 - p); nd[k + 1] += q * p
            dist = nd
        mu = sum(k * q for k, q in enumerate(dist))
        out.append((dist, mu, V))
    return out

def part2():
    print("\n== Part 2: Theorem D (single point) and Lemma E / Theorem F (per-group, per-bag) on G6a designs")
    for n, g, Xn, R in [(24, 4, 8, 2), (30, 5, 10, 2)]:
        F = random_3regular_formula(n, seed=n)
        X = list(range(1, Xn + 1)); M = 3 * n
        q = math.ceil((M * R + 1) / len(X)) - 1
        occ = occurrences(F)
        occ_of = {v: [o for o in occ if o[1] == v] for v in range(n)}
        for name, builder in [("cut-measure, variable-uniform", lambda: design_cut_measure(F, n, X, R, g, 1, True)),
                              ("cut-measure, clause-dependent", lambda: design_cut_measure(F, n, X, R, g, 2, False)),
                              ("global-sum", lambda: design_global_sum(F, n, X, R, g, 3)),
                              ("random dense", lambda: design_random_dense(F, n, X, R, g, 4))]:
            P, groups, alpha, Rsp = builder()
            st = check_design(F, n, X, R, P, groups, alpha, Rsp)
            rng = random.Random(n)
            omegas = [{v: rng.choice(occ_of[v]) for v in range(n)} for _ in range(3)]
            sumZ = 0; sumEN = 0.0; nul_tot = [0] * len(omegas); good = []; rows = []
            minslackE = float("inf"); betas = []
            for i, t in enumerate(X):
                gt = groups[t]; at = alpha[t]
                r = [rank_Q(Rsp[t][h]) for h in range(len(gt))]
                Zt = sum(1 for o in occ if Pval(P[o], t).get(o[1], 0) == 0)
                D = group_dists(F, n, gt, at)
                EN = 0.0; bagterm = {}
                for h, (dist, mu, V) in enumerate(D):
                    a = len(gt[h])
                    EDh = sum(qq * max(a - r[h] - k, 0) for k, qq in enumerate(dist))
                    EN += EDh
                    if mu <= a + 1e-12:   # Lemma E
                        lhs = EDh; rhs = C0 / 2 * min(V, math.sqrt(V)) - r[h]
                        assert lhs + 1e-9 >= rhs, ("Lemma E violated", t, h, lhs, rhs)
                        minslackE = min(minslackE, lhs - rhs)
                        bagterm[h] = C0 / 2 * min(V, math.sqrt(V))
                    if r[h] > 0:
                        betas.append((V / a, mu <= a + 1e-12))
                # strong Lemma 1 form: nullity >= sum_h (a_h - r_h - |I_h|)_+
                for k, om in enumerate(omegas):
                    cnt = [0] * len(gt)
                    for v in range(n):
                        cnt[at[om[v][0]]] += 1
                    Nk = sum(max(len(gt[h]) - r[h] - cnt[h], 0) for h in range(len(gt)))
                    nul = phi_matrix_nullity(F, n, t, P, om)
                    assert nul >= Nk, ("strong rank bound violated", t, nul, Nk)
                    nul_tot[k] += nul
                sumZ += Zt; sumEN += EN
                if Zt + EN < 2 * (q + 1):
                    good.append(t)
                # torso and reading-weighted width
                I = interaction_graph(F, t, P, gt, at)
                meas = {h for h in range(len(gt)) if r[h] > 0}
                Tor = torso(I, meas)
                W = weighted_tw(Tor.keys(), Tor, {h: r[h] for h in Tor})
                wt = weighted_tw(I.keys(), I, {h: 1 for h in I}) - 1
                Wa = weighted_tw(Tor.keys(), Tor, {h: len(gt[h]) for h in Tor})
                # Theorem F at this t (valid at every t with EN replacing 2q+2): every bag B of the
                # W-optimal decomposition has r(B) >= sum_{B,ns} c0/2 min(V,sqrtV) - EN. Check on the
                # whole measured set's best bag via brute force over cliques of a min-W triangulation is
                # heavy; check the weaker per-vertex-set form on every subset of measured groups:
                for k in range(1, len(meas) + 1):
                    for B in itertools.combinations(sorted(meas), k):
                        lhsB = sum(r[h] for h in B)
                        rhsB = sum(bagterm.get(h, 0.0) for h in B) - EN
                        assert lhsB + 1e-9 >= rhsB, ("Theorem F (subset form) violated", t, B)
                # Theorem C key step: every edge of Pr(F') (Z_t occurrences removed) joins equal or
                # I_t-adjacent groups, and no two unmeasured groups are adjacent (asserted in torso()).
                gid = {u: hh for hh, S in enumerate(gt) for u in S}
                for e, c in enumerate(F):
                    alive = [v for (v, _) in c if Pval(P[(e, v)], t).get(v, 0) != 0]
                    for u1, u2 in itertools.combinations(alive, 2):
                        assert gid[u1] == gid[u2] or gid[u2] in I[gid[u1]], ("Thm C blow-up fails", t, e)
                Gam = weighted_tw(Tor.keys(), Tor, {h: bagterm.get(h, 0.0) for h in Tor})
                assert W + 1e-9 >= Gam - EN, ("Theorem F violated", t, W, Gam, EN)
                # fraction of variables of measured groups that are non-surplus with V_h >= 0.15 a_h
                Dm = {h: D[h] for h in meas}
                good_a = sum(len(gt[h]) for h in meas if Dm[h][1] <= len(gt[h]) + 1e-12 and Dm[h][2] >= 0.15 * len(gt[h]))
                tot_a = sum(len(gt[h]) for h in meas)
                rows.append(f"t={t} |Z|={Zt} E N={EN:.2f} rho={sum(r)} w_t={wt} W_t={W} W^a_t={Wa} s={st[i]['s']}"
                            f"  Gamma_t(c0/2-weighted)={Gam:.2f}  W_t>=Gamma_t-EN ok; frac vars in (ns, V>=.15a) measured groups={good_a/tot_a:.2f}")
            for k in range(len(omegas)):
                assert nul_tot[k] <= n * R
            assert sumZ <= M * R
            assert sumZ + sumEN <= (M + n) * R + 1e-9
            assert good, "Theorem D: no good point"
            adj = primal_adj(F, n); lb = mmd_plus_lb(adj)
            # Theorem C at every t: tw(PrF) <= W^a + s - 1 + |Z_t|  -> check with lower bound lb
            print(f"  -- {name}: n={n} |X|={Xn} R={R} q={q}: sum|Z_t|={sumZ}<=MR={M*R}; sum E N_t={sumEN:.2f};"
                  f" sum(|Z|+EN)<= (M+n)R={(M+n)*R}; good points (|Z|+EN<2(q+1)): {good};"
                  f" min Lemma-E slack={minslackE:.3f}; nullity sums={nul_tot}<=nR={n*R}; mmd+ lb tw(PrF)={lb}")
            bm = [b for b, ns in betas if ns]
            print(f"     measured groups: non-surplus {len(bm)}/{len(betas)}, V_h/a_h over non-surplus measured:"
                  f" min={min(bm) if bm else float('nan'):.3f} mean={sum(bm)/max(1,len(bm)):.3f}")
            for rr in rows:
                print("    ", rr)

if __name__ == "__main__":
    part1()
    part2()
