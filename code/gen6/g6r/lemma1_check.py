# Independent check of G6a Lemma 1 / Theorem B step on random schemes built to satisfy Def 5.1(c)(3).
# P_o(T) = e_v + sum_{k=1..R} c_k T^k, constraints P_o(t)|_{S_h} in R_{t,h} for h != alpha_t(e) (linear in c).
import random, sympy as sp, itertools
T = sp.symbols('T')
def run(seed, n=9, m=7, R=2, npts=4, rdim=1):
    rng = random.Random(seed)
    while True:
        F = [rng.sample(range(n), 3) for _ in range(m)]
        if len({v for c in F for v in c}) == n: break
    occ = [(e, v) for e, c in enumerate(F) for v in c]
    X = list(range(1, npts+1))
    design = {}
    for t in X:
        perm = list(range(n)); rng.shuffle(perm)
        g = 3; groups = [perm[i::g] for i in range(g)]
        alpha = [rng.randrange(g) for _ in range(m)]
        Rsp = [sp.Matrix([[rng.randint(-2,2) for _ in S] for _ in range(min(rdim,len(S)))]) for S in groups]
        design[t] = (groups, alpha, Rsp)
    P = {}
    for (e, v) in occ:
        cs = sp.symbols(f'c0:{R*n}')
        vec = [ (1 if u==v else 0) + sum(cs[k*n+u]*T**(k+1) for k in range(R)) for u in range(n)]
        eqs = []
        for t in X:
            groups, alpha, Rsp = design[t]
            for h, S in enumerate(groups):
                if h == alpha[e]: continue
                # vec(t)|_S in rowspace(Rsp[h]): orthogonal to nullspace of Rsp[h]
                ns = Rsp[h].nullspace() if Rsp[h].rows else [sp.eye(len(S))[:,i] for i in range(len(S))]
                for z in ns:
                    eqs.append(sum(z[j]*vec[u].subs(T,t) for j,u in enumerate(S)))
        sol = sp.linsolve(eqs, cs)
        sol = list(sol)
        if not sol: return None
        sol = sol[0]; free = set().union(*[s.free_symbols for s in sol])
        sub = {f: rng.randint(-3,3) for f in free}
        vals = [s.subs(sub) for s in sol]
        P[(e,v)] = [sp.expand(x.subs(dict(zip(cs, vals)))) for x in vec]
    # check condition (3) and P(0)=e_v
    for (e,v),vec in P.items():
        assert [x.subs(T,0) for x in vec] == [1 if u==v else 0 for u in range(n)]
    occs_of = {v: [o for o in occ if o[1]==v] for v in range(n)}
    rho = {t: sum(design[t][2][h].rank() for h in range(3)) for t in X}
    worst = None
    for trial in range(6):
        om = {v: rng.choice(occs_of[v]) for v in range(n)}
        Phi = sp.Matrix([P[om[v]] for v in range(n)])
        det = sp.Poly(sp.expand(Phi.det()), T)
        assert det.eval(0) == 1 and det.degree() <= n*R
        tot = 0
        for t in X:
            groups, alpha, Rsp = design[t]
            null = n - Phi.subs(T,t).rank()
            I = [sum(1 for v in range(n) if alpha[om[v][0]]==h) for h in range(3)]
            dfc = sum(max(len(S)-I[h],0) for h,S in enumerate(groups))
            assert null >= dfc - rho[t], (null, dfc, rho[t])
            # multiplicity of root t
            mult = 0; p = det
            while p.eval(t) == 0: p = sp.Poly(sp.quo(p.as_expr(), T-t), T); mult += 1
            assert mult >= null
            tot += max(dfc - rho[t], 0)
        assert tot <= n*R
        worst = max(worst or 0, tot)
    return rho, worst
for s in range(6):
    print(s, run(s))
