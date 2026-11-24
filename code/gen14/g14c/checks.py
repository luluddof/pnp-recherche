"""G14c checks. Run: PYTHONDONTWRITEBYTECODE=1 python3 checks.py > checks.txt"""
import sys, random, itertools, math, time
sys.dont_write_bytecode = True
from lib import *
C0 = (1 - math.exp(-2)) / 4

def gfun(y): return min(y, math.sqrt(y)) if y > 0 else 0.0

def clause_cliques(F, A, groups, alpha):
    out = []
    for e, c in enumerate(F):
        cl = {alpha[e]}
        for (v, _) in c:
            row = A[(e, v)]
            for h, S in enumerate(groups):
                if h != alpha[e] and any(row[u] for u in S): cl.add(h)
        out.append(cl)
    return out

def tw_of(vertices, edges):
    vs = sorted(vertices); ix = {v: i for i, v in enumerate(vs)}
    adj = {i: set() for i in range(len(vs))}
    for a, b in edges:
        if a != b and a in ix and b in ix:
            adj[ix[a]].add(ix[b]); adj[ix[b]].add(ix[a])
    return treewidth_exact(len(vs), adj)

def design_params(F, n, A, groups, alpha, Rsp):
    """s_t, rho'_t, w_t = tw(I_t), tw(G_M')+1, row-treewidth+1 in the basis Rsp (upper bound on lw_t)."""
    g = len(groups)
    cl = clause_cliques(F, A, groups, alpha)
    Iedges = [(a, b) for c in cl for a in c for b in c if a < b]
    w = tw_of(range(g), Iedges)
    Mset = [h for h in range(g) if len(Rsp[h]) > 0]
    nb = {h: set() for h in range(g)}
    for a, b in Iedges: nb[a].add(b); nb[b].add(a)
    tor = [(a, b) for a, b in Iedges if a in Mset and b in Mset]
    for u in range(g):
        if u not in Mset:
            assert all(x in Mset for x in nb[u])          # unmeasured groups pairwise non-adjacent
            tor += [(a, b) for a in nb[u] for b in nb[u] if a < b]
    twM = tw_of(Mset, tor)
    V, adj, scope, reads = row_graph(F, n, A, groups, alpha, Rsp)
    assert len(V) <= 16
    rtw = treewidth_exact(len(V), adj) + 1
    return dict(s=max(len(S) for S in groups), rho=len(V), w=w, twM1=twM + 1, rtw=rtw)

# ------------------------------------------------------------------ interface-table evaluator
def phi(clause, vals):
    fal = Fr(1)
    for (v, pos), y in zip(clause, vals):
        fal *= (1 - y) if pos else y
    return 1 - fal

def random_tree(leaves, rng):
    nodes = [(l,) for l in leaves]; rng.shuffle(nodes)
    tree = {}
    while len(nodes) > 1:
        i, j = rng.sample(range(len(nodes)), 2)
        a, b = nodes[i], nodes[j]
        c = tuple(sorted(a + b)); tree[c] = (a, b)
        nodes = [x for k, x in enumerate(nodes) if k not in (i, j)] + [c]
    return nodes[0], tree

def interface_dp(F, n, A, root, tree):
    """Evaluate H(t) by the interface-table DP along the rooted binary tree. Leaves: ('c',e), ('v',u).
    Returns (H, list of (node, table_size, kappa))."""
    occ = occurrences(F)
    info = {}
    def setup(u):
        E = [l[1] for l in u if l[0] == 'c']; S = [l[1] for l in u if l[0] == 'v']
        rin = [o for o in occ if o[0] in E]; rout = [o for o in occ if o[0] not in E]
        info[u] = (E, S, rin, rout)
        if u in tree:
            for ch in tree[u]: setup(ch)
    setup(root)
    alph = {}
    def alphas(u):
        if u in alph: return alph[u]
        E, S, rin, rout = info[u]
        if u not in tree:
            if S:
                v = S[0]
                res = {tuple(Fr(0) for o in rout), tuple(A[o][v] for o in rout)}
            else:
                res = {tuple(Fr(0) for o in rout)}
        else:
            u1, u2 = tree[u]
            ro1, ro2 = info[u1][3], info[u2][3]
            i1 = {o: k for k, o in enumerate(ro1)}; i2 = {o: k for k, o in enumerate(ro2)}
            res = {tuple(a1[i1[o]] + a2[i2[o]] for o in rout) for a1 in alphas(u1) for a2 in alphas(u2)}
        alph[u] = res; return res
    memo = {u: {} for u in info}
    def T(u, beta):
        if beta in memo[u]: return memo[u][beta]
        E, S, rin, rout = info[u]
        if u not in tree:
            if S:
                res = {}
                v = S[0]
                for x in (0, 1):
                    a = tuple(A[o][v] * x for o in rout); res[a] = res.get(a, 0) + 1
            else:
                e = E[0]
                bi = {o: k for k, o in enumerate(rin)}
                res = {tuple(Fr(0) for o in rout): phi(F[e], [beta[bi[(e, v)]] for (v, _) in F[e]])}
        else:
            u1, u2 = tree[u]
            E1, S1, ri1, ro1 = info[u1]; E2, S2, ri2, ro2 = info[u2]
            bi = {o: k for k, o in enumerate(rin)}
            i1 = {o: k for k, o in enumerate(ro1)}; i2 = {o: k for k, o in enumerate(ro2)}
            res = {a: Fr(0) for a in alphas(u)}
            for a1 in alphas(u1):
                b2 = tuple(beta[bi[o]] + a1[i1[o]] for o in ri2)
                T2 = T(u2, b2)
                for a2 in alphas(u2):
                    b1 = tuple(beta[bi[o]] + a2[i2[o]] for o in ri1)
                    val = T(u1, b1)[a1] * T2[a2]
                    a = tuple(a1[i1[o]] + a2[i2[o]] for o in rout)
                    res[a] += val
        memo[u][beta] = res; return res
    H = T(root, tuple(Fr(0) for o in occ))[()]
    return H, memo, info

def brute_H(F, n, A):
    tot = Fr(0)
    for bits in range(1 << n):
        prod = Fr(1)
        for e, c in enumerate(F):
            vals = [sum((A[(e, v)][u] for u in range(n) if bits >> u & 1), Fr(0)) for (v, _) in c]
            prod *= phi(c, vals)
            if prod == 0: break
        tot += prod
    return tot

def make_designs(F, n, X, R, g, seed):
    out = {}
    out['cut_uniform'] = D.design_cut_measure(F, n, X, R, g, seed, True)
    out['cut_clausedep'] = D.design_cut_measure(F, n, X, R, g, seed, False)
    out['global_sum'] = D.design_global_sum(F, n, X, R, g, seed)
    out['random_dense'] = D.design_random_dense(F, n, X, R, g, seed, k=2)
    return out

def nullity(F, n, A, omega):
    return n - rank_Q([A[omega[v]] for v in range(n)])

# =================================================================== PART 1
def part1():
    print("PART 1  interface-table DP: correctness and table size >= 2^kappa  (soundness (c))")
    rng = random.Random(1); tot = 0; minratio = 1e9
    for n, seed in ((6, 1), (6, 2), (7, 3)):
        F = random_3regular_formula(n, seed); X = [1, 2, 3, 4, 5, 6]; R = 2
        for name, (P, groups, alpha, Rsp) in make_designs(F, n, X, R, 2, seed).items():
            for t in X[:2]:
                A = form_matrix(F, n, P, t); C = Cuts(F, n, A); Hb = brute_H(F, n, A)
                leaves = [('c', e) for e in range(len(F))] + [('v', u) for u in range(n)]
                for rep in range(2):
                    root, tree = random_tree(leaves, rng)
                    H, memo, info = interface_dp(F, n, A, root, tree)
                    assert H == Hb, (name, t, H, Hb)
                    for u in info:
                        if u == root: continue
                        E, S, rin, rout = info[u]
                        Em = sum(1 << e for e in E); Sm = sum(1 << v for v in S)
                        kap = C.kappa_ES(Em, Sm)
                        size = sum(len(d) for d in memo[u].values())
                        assert size >= 2 ** kap, (name, t, u, size, kap)
                        minratio = min(minratio, size / 2 ** kap); tot += 1
    print(f"  DP value == brute-force H(t) in all runs; {tot} node tables checked, size >= 2^kappa always"
          f" (min ratio size/2^kappa = {minratio:.2f})")

# =================================================================== PART 2
def part2():
    print("PART 2  cut inequality (*) and multi-cut inequality (Lemma 1/Thm 3), exact ranks")
    rng = random.Random(2)
    for n, seed in ((6, 11), (6, 12)):
        F = random_3regular_formula(n, seed); m = len(F); X = [1, 2, 3, 4, 5, 6]; R = 2
        occ_of = {v: [o for o in occurrences(F) if o[1] == v] for v in range(n)}
        omegas = [dict(zip(range(n), ch)) for ch in itertools.product(*[occ_of[v] for v in range(n)])]
        for name, (P, groups, alpha, Rsp) in make_designs(F, n, X, R, 2, seed).items():
            cnt = tight = 0; multi = 0; nulsum = [0] * len(omegas)
            for t in X:
                A = form_matrix(F, n, P, t); C = Cuts(F, n, A)
                nul = [nullity(F, n, A, om) for om in omegas]
                for k, x in enumerate(nul): nulsum[k] += x
                if t > 2: continue
                homes = [sum(1 << v for v in range(n) if False) for _ in omegas]
                # X_1(omega, Emask) = #{v : clause(omega(v)) in Emask}
                cl = [[om[v][0] for v in range(n)] for om in omegas]
                for Em in range(1 << m):
                    x1 = [sum(Em >> c[v] & 1 for v in range(n)) for c in cl]
                    for Sm in range(1 << n):
                        a = bin(Sm).count('1')
                        r12 = C.rk(Em, C.fullS ^ Sm); r21 = C.rk(C.fullE ^ Em, Sm)
                        for k in range(0, len(omegas), 7):
                            Dv = a - x1[k]
                            need = (Dv - r21) if Dv > 0 else (-Dv - r12)
                            assert nul[k] >= need, (name, t, Em, Sm, k)
                            cnt += 1; tight += nul[k] == need
                # multi-cut on random trees
                leaves = [('c', e) for e in range(m)] + [('v', u) for u in range(n)]
                for rep in range(30):
                    root, tree = random_tree(leaves, rng)
                    nodes = [u for u in tree for u in tree[u]]          # every non-root node = one tree edge
                    Y = rng.sample(nodes, rng.randint(1, 5))
                    # parts: components after cutting edges above nodes in Y
                    def part_of(leaf):
                        best = None
                        for u in Y:
                            if leaf in u and (best is None or len(u) < len(best)): best = u
                        return best
                    parts = {}
                    for l in leaves: parts.setdefault(part_of(l), []).append(l)
                    kapsum = 0
                    for u in Y:
                        Em = sum(1 << l[1] for l in u if l[0] == 'c'); Sm = sum(1 << l[1] for l in u if l[0] == 'v')
                        kapsum += C.kappa_ES(Em, Sm)
                    for k in range(0, len(omegas), 13):
                        om = omegas[k]; dfc = 0; strong = 0
                        for key, ls in parts.items():
                            Es = {l[1] for l in ls if l[0] == 'c'}; Ss = [l[1] for l in ls if l[0] == 'v']
                            xi = sum(1 for v in range(n) if om[v][0] in Es)
                            dfc += max(len(Ss) - xi, 0)
                            Em = sum(1 << e for e in Es); Sm = sum(1 << v for v in Ss)
                            strong += max(len(Ss) - xi - C.rk(C.fullE ^ Em, Sm), 0)
                        assert nul[k] >= strong, "partition rank lemma"
                        assert kapsum + nul[k] >= dfc, "multi-cut"
                        multi += 1
            assert max(nulsum) <= n * R
            print(f"  n={n} {name:14s}: (*) checked on {cnt} (cut,omega) pairs at t=1,2 (all 4096 cuts), tight in {tight};"
                  f" multi-cut: {multi} cases; max_omega sum_t nullity = {max(nulsum)} <= nR = {n*R} (all {len(omegas)} omegas)")

# =================================================================== PART 3
def part3():
    print("PART 3  exact branch-width bw_t (IR model) vs DL parameters and the row graph")
    print("   columns: s_t, rho'_t, w_t=tw(I_t), tw(G_M')+1, rtw=row-treewidth+1 (natural reduced basis), bw_t (exact)")
    for n, seed, g in ((6, 21, 2), (6, 22, 3), (7, 23, 2), (7, 24, 3)):
        F = random_3regular_formula(n, seed); X = [1, 2, 3, 4, 5, 6]; R = 2
        twP = treewidth_exact(n, primal_adj(F, n))
        # trivial scheme (paper Prop. 6.1): P_o = e_v, fully measured consecutive blocks
        A = form_matrix(F, n, trivial_P(F, n), 1); C = Cuts(F, n, A); t0 = time.time(); bw0 = branchwidth(C)
        s = math.ceil(math.sqrt(n)); gs = [list(range(i, min(i + s, n))) for i in range(0, n, s)]
        assert 2 * (twP + 1) <= 3 * bw0, "non-degeneracy (b) violated"
        print(f"  n={n} seed={seed}: tw(Pr F)={twP};  TRIVIAL scheme: bw={bw0}  [(b): (2/3)(tw+1)={2*(twP+1)/3:.2f} <= bw],"
              f"  DL DP-type exponent = ceil(sqrt n) = {s}   ({time.time()-t0:.0f}s)")
        for name, (P, groups, alpha, Rsp) in make_designs(F, n, X, R, g, seed).items():
            for t in X[:2]:
                A = form_matrix(F, n, P, t); C = Cuts(F, n, A)
                p = design_params(F, n, A, groups[t], alpha[t], Rsp[t])
                bw = branchwidth(C)
                assert p['w'] <= p['twM1'] <= p['rtw'], ("(d) chain violated", p)
                assert bw <= p['rtw'] + p['s'], ("bw <= lw + s violated", bw, p)
                assert p['rtw'] <= p['rho']
                print(f"     {name:14s} t={t}: s={p['s']} rho'={p['rho']} w={p['w']} tw(G_M')+1={p['twM1']} rtw={p['rtw']} bw={bw}")

# =================================================================== PART 4
def part4():
    print("PART 4  aggregate-reading example: the row model (A) is NOT sound for evaluators with derived variables")
    for mm, s in ((4, 2), (6, 2), (8, 3), (12, 3)):
        # groups h_1..h_mm of s variables; hub variable u; hub clause reads sum_i 1_{S_i}; each group has one internal clause
        n = mm * s + 1; u = n - 1
        F = [[(u, True)]] + [[(i * s + j, True) for j in range(s)] for i in range(mm)]
        A = {(0, u): [Fr(1)] * (n - 1) + [Fr(1)]}
        for i in range(mm):
            for j in range(s):
                row = [Fr(0)] * n; row[i * s + j] = Fr(1); A[(i + 1, i * s + j)] = row
        C = Cuts(F, n, A)
        # row model: hub clique contains one functional per satellite group: lw >= mm (basis-free: sum_h dim L_{g->h})
        lw_lb = mm
        # explicit caterpillar: group subtree (clause i+1, its variables), chained; hub leaves at the end
        widths = []
        acc_E = 0; acc_S = 0
        for i in range(mm):
            gE = 1 << (i + 1); gS = 0
            widths.append(C.kappa_ES(gE, 0))
            for j in range(s):
                widths.append(C.kappa_ES(0, 1 << (i * s + j)))
                gS |= 1 << (i * s + j); widths.append(C.kappa_ES(gE, gS))
            acc_E |= gE; acc_S |= gS; widths.append(C.kappa_ES(acc_E, acc_S))
        widths += [C.kappa_ES(1, 0), C.kappa_ES(0, 1 << u), C.kappa_ES(1, 1 << u)]
        print(f"  {mm} satellite groups x {s} vars: row-model width lw >= {lw_lb} (charge 2^{lw_lb});"
              f" IR model: explicit tree of width {max(widths)}; honest partial-sum table has <= {mm*s+1} values per cut")

# =================================================================== PART 5
def part5():
    print("PART 5  Theorem 2 ingredients: balanced edge / balanced node always exists; anticoncentration on cuts")
    rng = random.Random(5); case1 = case2 = 0; worst = 1e9
    for n, seed in ((6, 31), (7, 32), (9, 33), (12, 34)):
        F = random_3regular_formula(n, seed); m = len(F)
        occ_cl = {v: [e for e, c in enumerate(F) for (x, _) in c if x == v] for v in range(n)}
        leaves = [('c', e) for e in range(m)] + [('v', u) for u in range(n)]
        for rep in range(200):
            root, tree = random_tree(leaves, rng)
            def J(Es): return sum(1 for v in range(n) if all(e in Es for e in occ_cl[v]))
            allE = set(range(m)); found = None
            for u in tree:
                for ch in tree[u]:
                    Es = {l[1] for l in ch if l[0] == 'c'}
                    if 3 * J(Es) <= 2 * n and 3 * J(allE - Es) <= 2 * n: found = ('edge', Es)
            if found:
                case1 += 1; Es = found[1]
                ps = [sum(e in Es for e in occ_cl[v]) / len(occ_cl[v]) for v in range(n)]
                dist = [1.0]
                for p in ps:
                    nd = [0.0] * (len(dist) + 1)
                    for k, q in enumerate(dist): nd[k] += q * (1 - p); nd[k + 1] += q * p
                    dist = nd
                V = sum(p * (1 - p) for p in ps)
                for a in range(n + 1):
                    ea = sum(q * abs(k - a) for k, q in enumerate(dist))
                    assert ea >= C0 * gfun(V) - 1e-12; worst = min(worst, ea - C0 * gfun(V))
            else:
                ok = False
                for u in tree:        # internal node: branches = two children + complement
                    ch = tree[u]
                    br = [{l[1] for l in c if l[0] == 'c'} for c in ch] + [allE - {l[1] for l in u if l[0] == 'c'}]
                    if all(3 * J(b) < n for b in br): ok = True
                assert ok, "neither a balanced edge nor a balanced node"
                case2 += 1
    print(f"  800 random trees: balanced edge (case 1) in {case1}, balanced node (case 2) in {case2}, none: 0;"
          f" anticoncentration E|X-a| >= c0*min(V,sqrt V) on all case-1 cuts (min slack {worst:.3f})")
    # value-count lemma: a rank-r matrix takes >= 2^r values on the cube
    bad = 0
    for it in range(300):
        s = rng.randint(1, 7); mrow = rng.randint(1, 5)
        B = [[Fr(rng.choice([0, 0, 1, -1, 2])) for _ in range(s)] for _ in range(mrow)]
        vals = {tuple(sum(b[i] for i in range(s) if x >> i & 1) for b in B) for x in range(1 << s)}
        assert len(vals) >= 2 ** rank_Q(B)
    print("  value-count lemma |B{0,1}^s| >= 2^rank(B): 300 random integer matrices pass")


# =================================================================== PART 6
def part6():
    print("PART 6  DL's 2^{w_t} is design-dependent: path formula, trivial fully-measured scheme, scattered groups")
    for n, lab in ((8, [0, 1, 2, 3, 1, 3, 0, 2]), (15, [0, 1, 2, 3, 4, 0, 2, 4, 1, 3, 0, 3, 1, 4, 2])):
        F = [[(i, True), (i + 1, False)] for i in range(n - 1)]
        A = form_matrix(F, n, trivial_P(F, n), 1); C = Cuts(F, n, A)
        g = max(lab) + 1; groups = [[v for v in range(n) if lab[v] == h] for h in range(g)]
        alpha = [lab[c[0][0]] for c in F]
        Rsp = [[[Fr(int(i == j)) for i in range(len(S))] for j in range(len(S))] for S in groups]
        cl = clause_cliques(F, A, groups, alpha)
        w = tw_of(range(g), [(a, b) for c in cl for a in c for b in c if a < b])
        widths = []; accE = accS = 0
        for i in range(n):
            accS |= 1 << i; widths += [C.kappa_ES(0, 1 << i), C.kappa_ES(accE, accS)]
            if i < n - 1:
                accE |= 1 << i; widths += [C.kappa_ES(1 << i, 0), C.kappa_ES(accE, accS)]
        twP = treewidth_exact(n, primal_adj(F, n))
        print(f"  path n={n}, {g} groups: w_t = tw(I_t) = {w}, s_t = {max(map(len, groups))};  IR: explicit caterpillar width {max(widths)}"
              f" (tw(Pr F) = {twP}); so 2^(w_t) is not a lower bound for DP evaluators here, 2^(bw_t) <= 2^{max(widths)}")

if __name__ == '__main__':
    t0 = time.time()
    for f in (part1, part2, part3, part4, part5, part6):
        f(); print(f"   [{time.time()-t0:.0f}s]"); sys.stdout.flush()
    print("ALL CHECKS PASSED")
