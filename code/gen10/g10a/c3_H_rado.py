# Lemma H (own construction, exact rank mod p) and tightness of nullity >= N^omega for generic rows.
import random, itertools
P = 2**31 - 1
def rank_mod(M):
    M = [r[:] for r in M]; rk = 0; ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(M)) if M[i][c] % P), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]
        inv = pow(M[rk][c], P-2, P)
        M[rk] = [x*inv % P for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c] % P:
                f = M[i][c]; M[i] = [(x - f*y) % P for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk
def formula(n, rng):
    while True:
        pool = [v for v in range(n) for _ in range(3)]; rng.shuffle(pool)
        cl = [tuple(pool[3*i:3*i+3]) for i in range(n)]
        if all(len(set(e)) == 3 for e in cl): return cl
def point(n, g, rng, dense):
    cl = formula(n, rng)
    lab = [rng.randrange(g) for _ in range(n)]
    S = [[v for v in range(n) if lab[v] == h] for h in range(g)]
    r = [rng.randint(0, min(2, len(S[h]))) for h in range(g)]
    # basis of R_h: r_h random vectors supported on S_h
    Rb = [[[ (rng.randrange(1, P) if v in S[h] else 0) for v in range(n)] for _ in range(r[h])] for h in range(g)]
    alpha = [lab[rng.choice(e)] if rng.random() < .8 else rng.randrange(g) for e in cl]
    rows = {}; adj = {h: set() for h in range(g)}
    for e, vs in enumerate(cl):
        clique = {alpha[e]}
        for v in vs:
            vec = [0]*n
            for u in S[alpha[e]]: vec[u] = rng.randrange(P) if rng.random() < .7 else 0
            for h in range(g):
                if h == alpha[e] or r[h] == 0: continue
                if dense or lab[v] == h or rng.random() < .15:
                    co = [rng.randrange(P) for _ in range(r[h])]
                    if rng.random() < .1: co = [0]*r[h]
                    add = [sum(c*b[u] for c, b in zip(co, Rb[h])) % P for u in range(n)]
                    if any(add): clique.add(h)
                    vec = [(x+y) % P for x, y in zip(vec, add)]
            rows[(e, v)] = vec
        for x in clique:
            for y in clique:
                if x != y: adj[x].add(y)
    return cl, lab, S, r, alpha, rows, adj
rng = random.Random(7)
nchk = tight = 0; rado_eq = rado_tot = 0
for it in range(120):
    n = rng.choice([9, 12, 15]); g = rng.choice([3, 4, 5]); dense = it % 2 == 0
    cl, lab, S, r, alpha, rows, adj = point(n, g, rng, dense)
    occ = {v: [e for e, vs in enumerate(cl) if v in vs] for v in range(n)}
    a = [len(x) for x in S]
    for _ in range(3):
        om = {v: rng.choice(occ[v]) for v in range(n)}
        Phi = [rows[(om[v], v)] for v in range(n)]
        nul = n - rank_mod(Phi)
        X = [sum(1 for v in range(n) if alpha[om[v]] == h) for h in range(g)]
        N = sum(max(a[h]-r[h]-X[h], 0) for h in range(g))
        assert nul >= N
        if dense:
            rado_tot += 1; rado_eq += (nul == N)
        for k in range(1, g+1):
            for A in itertools.combinations(range(g), k):
                A = set(A); dA = set().union(*(adj[h] for h in A)) - A
                bound = sum(X[h]-a[h] for h in A) - sum(r[h] for h in dA)
                assert nul >= bound, (nul, bound)
                nchk += 1; tight += (nul == bound and bound > 0)
print(f"Lemma H: {nchk} (point,omega,A) checks pass, {tight} tight with positive bound")
print(f"dense interaction, 'generic' rows: nullity == N^omega in {rado_eq}/{rado_tot} cases (note: alpha-block entries 70% dense, coefficients 10% zero)")
