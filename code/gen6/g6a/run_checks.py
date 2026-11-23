"""Exact checks of designs + numerical check of Theorems A/B inequality chain."""
import sys, random, math, itertools, json
from fractions import Fraction as Fr
from scheme import *
from designs import *

C0 = 0.216

def analyse(name, F, n, X, R, P, groups, alpha, Rsp, n_omega=4, seed=0):
    M = 3 * n
    st = check_design(F, n, X, R, P, groups, alpha, Rsp)   # raises if (3) fails
    ws = interaction_width(F, X, P, groups, alpha)
    q = math.ceil((M * R + 1) / len(X)) - 1                 # uniform Taylor order, sum(q+1) >= MR+1
    rng = random.Random(seed)
    occ_of = {v: [o for o in occurrences(F) if o[1] == v] for v in range(n)}
    lhs_sum = 0.0
    rows = []
    worst_gap = None
    nul_tot = {k: 0 for k in range(n_omega)}
    omegas = [{v: rng.choice(occ_of[v]) for v in range(n)} for _ in range(n_omega)]
    for i, t in enumerate(X):
        Ed, split, imb = expected_deficiency(F, n, groups[t], alpha[t])
        s = st[i]['s']; rho = st[i]['rho']
        lb = max(C0 * 2 / (9 * math.sqrt(2)) * split / math.sqrt(s), imb)
        assert Ed + 1e-9 >= lb, ("Theorem B lower bound on E def violated", Ed, lb)
        for k, om in enumerate(omegas):
            nul = phi_matrix_nullity(F, n, t, P, om)
            d = deficiency(groups[t], alpha[t], om, n)
            assert nul >= d - rho, ("rank lemma violated", t, nul, d, rho)
            nul_tot[k] += nul
        lhs_sum += max(Ed - rho, 0)
        rows.append(dict(t=t, s=s, g=st[i]['g'], rho=rho, mu=st[i]['mu'], w_ub=ws[i], q=q,
                         split=split, Edef=round(Ed, 2), lb=round(lb, 2)))
    for k in nul_tot:
        assert nul_tot[k] <= n * R, ("sum of nullities exceeds deg det", nul_tot[k], n * R)
    assert lhs_sum <= n * R + 1e-9
    print(f"== {name}: n={n} |X|={len(X)} R={R} q={q}  sum_t(Edef-rho)_+={lhs_sum:.2f} <= nR={n*R}"
          f"  max nullity-sum over omegas={max(nul_tot.values())}")
    for r in rows:
        print("   ", r)
    return rows

# ------------------------------------------------------------- end-to-end exact evaluation
def brute_H_at(F, n, P, t):
    vals = {o: Pval(P[o], t) for o in P}
    tot = Fr(0)
    for bits in range(1 << n):
        x = [(bits >> u) & 1 for u in range(n)]
        prod = Fr(1)
        for e, c in enumerate(F):
            fal = Fr(1)
            for (v, pos) in c:
                l = sum((a for u, a in vals[(e, v)].items() if x[u]), Fr(0))
                fal *= (1 - l) if pos else l
            prod *= (1 - fal)
            if prod == 0:
                break
        tot += prod
    return tot

def design_H_at(F, n, P, groups_t, alpha_t, Rsp_t, t):
    """Product-type evaluation: enumerate reading tuples, count inside each group."""
    g = len(groups_t)
    B = [Rsp_t[h] for h in range(g)]
    # readings of each group
    def reading(h, xs):
        return tuple(sum((b[i] for i, xi in enumerate(xs) if xi), Fr(0)) for b in B[h])
    Yset = []
    for h, S in enumerate(groups_t):
        Yset.append(sorted({reading(h, xs) for xs in itertools.product((0, 1), repeat=len(S))}))
    # coefficients lambda: foreign restriction = sum_j lam_j * B[h][j]
    def coords(h, vec):
        # solve vec = sum lam_j B[h][j] (B rows independent, reduced) by least elimination
        rows = B[h]; k = len(rows)
        import sympy
        A = sympy.Matrix([[rows[j][i] for j in range(k)] for i in range(len(vec))])
        sol = A.solve_least_squares(sympy.Matrix(vec)) if k else sympy.Matrix([])
        lam = [Fr(int(sympy.fraction(s)[0]), int(sympy.fraction(s)[1])) for s in sol]
        assert all(sum((lam[j] * rows[j][i] for j in range(k)), Fr(0)) == vec[i] for i in range(len(vec)))
        return lam
    lin = {}
    for e, c in enumerate(F):
        ge = alpha_t[e]
        for (v, pos) in c:
            val = Pval(P[(e, v)], t)
            own = [val.get(u, Fr(0)) for u in groups_t[ge]]
            fo = {}
            for h, S in enumerate(groups_t):
                if h == ge:
                    continue
                vec = [val.get(u, Fr(0)) for u in S]
                if any(vec):
                    fo[h] = coords(h, vec)
            lin[(e, v)] = (own, fo)
    clauses_of = [[e for e in range(len(F)) if alpha_t[e] == h] for h in range(g)]
    total = Fr(0); ntuples = 0
    for ys in itertools.product(*Yset):
        ntuples += 1
        prod = Fr(1)
        for h, S in enumerate(groups_t):
            part = Fr(0)
            for xs in itertools.product((0, 1), repeat=len(S)):
                if reading(h, xs) != ys[h]:
                    continue
                pr = Fr(1)
                for e in clauses_of[h]:
                    fal = Fr(1)
                    for (v, pos) in F[e]:
                        own, fo = lin[(e, v)]
                        l = sum((a for a, xi in zip(own, xs) if xi), Fr(0))
                        for hh, lam in fo.items():
                            l += sum((lj * yj for lj, yj in zip(lam, ys[hh])), Fr(0))
                        fal *= (1 - l) if pos else l
                    pr *= (1 - fal)
                part += pr
            prod *= part
            if prod == 0:
                break
        total += prod
    return total, ntuples, [len(Y) for Y in Yset]

PR = (1 << 61) - 1
def H_poly_modp(F, n, P, R):
    M = 3 * n
    def mp(c):
        return [(x.numerator * pow(x.denominator, PR - 2, PR)) % PR for x in c]
    Pm = {o: {u: mp(c) for u, c in P[o].items()} for o in P}
    def mul(a, b):
        r = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    r[i + j] = (r[i + j] + x * y) % PR
        return r
    H = [0] * (M * R + 1)
    for bits in range(1 << n):
        x = [(bits >> u) & 1 for u in range(n)]
        prod = [1]
        for e, c in enumerate(F):
            fal = [1]
            for (v, pos) in c:
                l = [0] * (R + 1)
                for u, cc in Pm[(e, v)].items():
                    if x[u]:
                        for i, a in enumerate(cc):
                            l[i] = (l[i] + a) % PR
                lit = [(-a) % PR for a in l]; lit[0] = (lit[0] + 1) % PR
                fal = mul(fal, lit if pos else l)
            phi = [(-a) % PR for a in fal]; phi[0] = (phi[0] + 1) % PR
            prod = mul(prod, phi)
        for i, a in enumerate(prod):
            H[i] = (H[i] + a) % PR
    return H

def hermite_recover(H, X, q, D):
    """Taylor data of H at X up to order q, then solve for deg<=D poly and evaluate at 0 (mod p)."""
    eqs = []
    for t in X:
        for j in range(q + 1):
            # j-th Taylor coefficient functional at t: sum_k C(k,j) t^(k-j) a_k
            row = [(math.comb(k, j) * pow(t % PR, k - j, PR)) % PR if k >= j else 0 for k in range(D + 1)]
            rhs = sum(r * a for r, a in zip(row, H + [0] * (D + 1 - len(H)))) % PR
            eqs.append(row + [rhs])
    # Gaussian elimination mod p
    m = len(eqs); col = 0; r = 0; piv_cols = []
    A = eqs
    for col in range(D + 1):
        piv = next((i for i in range(r, m) if A[i][col] % PR), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = pow(A[r][col], PR - 2, PR)
        A[r] = [(a * inv) % PR for a in A[r]]
        for i in range(m):
            if i != r and A[i][col]:
                f = A[i][col]
                A[i] = [(a - f * b) % PR for a, b in zip(A[i], A[r])]
        piv_cols.append(col); r += 1
    assert len(piv_cols) == D + 1, "Hermite data do not determine degree-D polynomials"
    return A[piv_cols.index(0)][-1]

def count_sat(F, n):
    c = 0
    for bits in range(1 << n):
        if all(any(((bits >> v) & 1) == pos for (v, pos) in cl) for cl in F):
            c += 1
    return c

if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "theorem":
        for n, g, Xn, R in [(24, 4, 8, 2), (30, 5, 10, 2), (36, 6, 9, 3)]:
            F = random_3regular_formula(n, seed=n)
            X = list(range(1, Xn + 1))
            for name, builder in [("cut-measure, variable-uniform", lambda: design_cut_measure(F, n, X, R, g, 1, True)),
                                  ("cut-measure, clause-dependent", lambda: design_cut_measure(F, n, X, R, g, 2, False)),
                                  ("global-sum (1 all-ones row/group + leftovers)", lambda: design_global_sum(F, n, X, R, g, 3)),
                                  ("random dense deformation, minimal R_{t,h}", lambda: design_random_dense(F, n, X, R, g, 4))]:
                P, groups, alpha, Rsp = builder()
                analyse(name, F, n, X, R, P, groups, alpha, Rsp, seed=n)
    elif mode == "e2e":
        n, g, R = 12, 3, 1
        X = [1, 2, 3, 4, 5, 6]
        F = random_3regular_formula(n, seed=7)
        nsat = count_sat(F, n)
        for name, builder in [("cut-measure clause-dependent", lambda: design_cut_measure(F, n, X, R, g, 2, False)),
                              ("global-sum", lambda: design_global_sum(F, n, X, R, g, 3))]:
            P, groups, alpha, Rsp = builder()
            st = check_design(F, n, X, R, P, groups, alpha, Rsp)
            H = H_poly_modp(F, n, P, R)
            assert H[0] == nsat
            q = math.ceil((3 * n * R + 1) / len(X)) - 1
            rec = hermite_recover(H, X, q, 3 * n * R)
            print(f"== e2e {name}: #F={nsat}, H(0) mod p={H[0]}, Hermite-recovered={rec}, q={q}")
            for i, t in enumerate(X):
                bh = brute_H_at(F, n, P, t)
                dh, ntup, ysz = design_H_at(F, n, P, groups[t], alpha[t], Rsp[t], t)
                hp = sum(a * pow(t, k, PR) for k, a in enumerate(H)) % PR
                bmod = (bh.numerator * pow(bh.denominator, PR - 2, PR)) % PR
                assert bh == dh and hp == bmod, (t, bh, dh)
                print(f"   t={t}: H(t)={bh} (brute) == design evaluation; s={st[i]['s']} rho={st[i]['rho']}"
                      f" reading-tuples={ntup} per-group readings={ysz}")
