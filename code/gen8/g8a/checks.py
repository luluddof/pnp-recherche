"""G8a checks. Imports G6a's scheme.py/designs.py read-only."""
import sys, random, math, itertools
sys.path.insert(0, '/home/user/pnp-recherche/code/gen6/g6a')
from fractions import Fraction as Fr
from scheme import *
from designs import *

# ---------- 1. anticoncentration used in Theorem B2: E|X-X'| >= 0.378 sqrt(N), X,X' iid Bin(N,1/2)
def EabsZ(N):
    # Z = Bin(2N,1/2) - N ; exact with integers
    tot = sum(math.comb(2*N, k) * abs(k - N) for k in range(2*N + 1))
    return Fr(tot, 2**(2*N))
worst = min((float(EabsZ(N)) / math.sqrt(N), N) for N in range(1, 1201))
print("1. min_N<=1200 E|X-X'|/sqrt(N) =", round(worst[0], 5), "at N =", worst[1], "(claim >= 0.378)")
assert worst[0] >= 0.378

# ---------- 2. Lemma H exact: nullity Phi^w(t) >= sum_{h in A}(X_h - a_h) - r(dA)
def interaction_adj(F, P, groups_t, alpha_t, t):
    gid = {u: h for h, S in enumerate(groups_t) for u in S}
    adj = {h: set() for h in range(len(groups_t))}
    for e, c in enumerate(F):
        cl = {alpha_t[e]}
        for (v, _) in c:
            for u in Pval(P[(e, v)], t):
                if gid[u] != alpha_t[e]:
                    cl.add(gid[u])
        for a, b in itertools.combinations(cl, 2):
            adj[a].add(b); adj[b].add(a)
    return adj

def lemmaH(n, seed, Xn, R, g, nom=3):
    F = random_3regular_formula(n, seed); X = list(range(1, Xn + 1))
    fams = [("cut-measure uniform", lambda: design_cut_measure(F, n, X, R, g, 1, True)),
            ("cut-measure clause-dep", lambda: design_cut_measure(F, n, X, R, g, 2, False)),
            ("global-sum", lambda: design_global_sum(F, n, X, R, g, 3)),
            ("random dense", lambda: design_random_dense(F, n, X, R, g, 4))]
    rng = random.Random(seed)
    occ_of = {v: [o for o in occurrences(F) if o[1] == v] for v in range(n)}
    for name, b in fams:
        P, groups, alpha, Rsp = b()
        check_design(F, n, X, R, P, groups, alpha, Rsp)
        nchk = 0; minslack = 10**9; tight = 0
        for t in X:
            G = len(groups[t]); a = [len(S) for S in groups[t]]
            r = [rank_Q(Rsp[t][h]) for h in range(G)]
            adj = interaction_adj(F, P, groups[t], alpha[t], t)
            for _ in range(nom):
                om = {v: rng.choice(occ_of[v]) for v in range(n)}
                nul = phi_matrix_nullity(F, n, t, P, om)
                Xh = [0] * G
                for v in range(n): Xh[alpha[t][om[v][0]]] += 1
                for k in range(1, G + 1):
                    for A in itertools.combinations(range(G), k):
                        Aset = set(A); dA = set().union(*(adj[h] for h in A)) - Aset
                        lb = sum(Xh[h] - a[h] for h in A) - sum(r[h] for h in dA)
                        assert nul >= lb, ("Lemma H violated", name, t, A, nul, lb)
                        nchk += 1; minslack = min(minslack, nul - lb); tight += (nul == lb)
        print(f"2. Lemma H  n={n} {name:24s}: {nchk} (t,omega,A) checks pass, min slack {minslack}, tight {tight}")

lemmaH(24, 7, 8, 2, 4)
lemmaH(30, 11, 10, 2, 5)

# ---------- 3. Shadow (Prop. S): iid labels, alpha(e)=label of a random var of e, fixed omega
def shadow(n, s, trials, seed):
    rng = random.Random(seed)
    F = random_3regular_formula(n, seed)
    b = 3; g = math.ceil(n / s); L = math.log(n)
    r = math.ceil(9 * math.sqrt(b * s * L))
    occ_of = {v: [e for e, c in enumerate(F) for (u, _) in c if u == v] for v in range(n)}
    worst = 0; maxvar = 0; Nsum = 0; Edef = 0; smax = 0
    for _ in range(trials):
        om = {v: rng.choice(occ_of[v]) for v in range(n)}
        lab = [rng.randrange(g) for _ in range(n)]
        al = [lab[rng.choice(c)[0]] for c in F]
        a = [0] * g; Xh = [0] * g
        for v in range(n): a[lab[v]] += 1; Xh[al[om[v]]] += 1
        dev = [a[h] - Xh[h] for h in range(g)]
        worst = max(worst, max(dev) / math.sqrt(b * s * L))
        Nsum += sum(max(a[h] - r - Xh[h], 0) for h in range(g))
        Edef += sum(max(d, 0) for d in dev) / trials
        smax = max(smax, max(a))
    print(f"3. shadow n={n} s={s} g={g}: r_h={r}, rho={g*r}, max s_t={smax}, mean def={Edef:.1f} "
          f"(bound sqrt(6bs)/2*g={math.sqrt(6*b*s)/2*g:.1f}), max (a-X)/sqrt(bs ln n)={worst:.3f} (<9), sum N^w={Nsum}")
    assert Nsum == 0
shadow(3000, round(3000 ** (2 / 3)), 200, 1)
shadow(8000, round(8000 ** (2 / 3)), 60, 2)
