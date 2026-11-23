import random
from fractions import Fraction as Fr
from scheme import occurrences

def polymul(a, b):
    r = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return r

def cutpoly(taus):
    c = [Fr(1)]
    for tau in taus:
        c = polymul(c, [Fr(1), Fr(-1, tau)])  # (1 - T/tau)
    return c

def random_partitions(n, X, g, rng):
    groups = {}
    for t in X:
        perm = list(range(n)); rng.shuffle(perm)
        groups[t] = [sorted(perm[h::g]) for h in range(g)]
    return groups

def majority_alpha(F, X, groups, rng):
    alpha = {}
    for t in X:
        gid = {u: h for h, S in enumerate(groups[t]) for u in S}
        a = []
        for c in F:
            cnt = {}
            for (v, _) in c:
                cnt[gid[v]] = cnt.get(gid[v], 0) + 1
            m = max(cnt.values())
            a.append(rng.choice([h for h, k in cnt.items() if k == m]))
        alpha[t] = a
    return alpha

def foreign(t, o, groups, alpha):
    e, v = o
    gid = {u: h for h, S in enumerate(groups[t]) for u in S}
    return gid[v] != alpha[t][e]

def minimal_measurements(F, n, X, P, groups, alpha):
    """R_{t,h} := span of all required restrictions (the smallest legal choice)."""
    from scheme import Pval
    Rsp = {}
    for t in X:
        Rsp[t] = []
        for h, S in enumerate(groups[t]):
            vecs = []
            for o in occurrences(F):
                if alpha[t][o[0]] == h:
                    continue
                val = Pval(P[o], t)
                vec = [val.get(u, Fr(0)) for u in S]
                if any(x != 0 for x in vec):
                    vecs.append(vec)
            Rsp[t].append(reduce_basis(vecs))
    return Rsp

def reduce_basis(vecs):
    M = [list(r) for r in vecs]
    basis = []
    for r in M:
        for b, piv in basis:
            if r[piv] != 0:
                f = r[piv] / b[piv]
                r = [x - f * y for x, y in zip(r, b)]
        nz = next((i for i, x in enumerate(r) if x != 0), None)
        if nz is not None:
            basis.append((r, nz))
    return [b for b, _ in basis]

def design_cut_measure(F, n, X, R, g, seed, uniform):
    """Each variable (uniform) or occurrence (clause-dependent) is cut at R points; measure the rest."""
    rng = random.Random(seed)
    groups = random_partitions(n, X, g, rng)
    alpha = majority_alpha(F, X, groups, rng)
    occ = occurrences(F)
    P = {}
    if uniform:
        for v in range(n):
            # choose R points where v has the most foreign occurrences
            score = [(sum(foreign(t, o, groups, alpha) for o in occ if o[1] == v), rng.random(), t) for t in X]
            taus = [t for _, _, t in sorted(score, reverse=True)[:R]]
            c = cutpoly(taus)
            for o in occ:
                if o[1] == v:
                    P[o] = {v: c}
    else:
        for o in occ:
            fpts = [t for t in X if foreign(t, o, groups, alpha)]
            rng.shuffle(fpts)
            taus = (fpts + [t for t in X if t not in fpts])[:R]
            P[o] = {o[1]: cutpoly(taus)}
    Rsp = minimal_measurements(F, n, X, P, groups, alpha)
    return P, groups, alpha, Rsp

def design_global_sum(F, n, X, R, g, seed):
    """P_o = d_o(T) e_v + gamma_o T * 1 ; all-ones row in every group + e_v where d_o(t) != 0."""
    rng = random.Random(seed)
    groups = random_partitions(n, X, g, rng)
    alpha = majority_alpha(F, X, groups, rng)
    P = {}
    for o in occurrences(F):
        fpts = [t for t in X if foreign(t, o, groups, alpha)]
        rng.shuffle(fpts)
        taus = (fpts + [t for t in X if t not in fpts])[:R]
        d = cutpoly(taus)
        gam = Fr(rng.randint(1, 5))
        Po = {u: [Fr(0), gam] for u in range(n)}
        Po[o[1]] = [d[0], d[1] + gam] + d[2:]
        P[o] = Po
    Rsp = minimal_measurements(F, n, X, P, groups, alpha)
    return P, groups, alpha, Rsp

def design_random_dense(F, n, X, R, g, seed, k=3):
    rng = random.Random(seed)
    groups = random_partitions(n, X, g, rng)
    alpha = majority_alpha(F, X, groups, rng)
    P = {}
    for o in occurrences(F):
        Po = {o[1]: [Fr(1)] + [Fr(0)] * R}
        for i in range(1, R + 1):
            for u in rng.sample(range(n), k):
                Po.setdefault(u, [Fr(0)] * (R + 1))
                Po[u][i] += rng.randint(-3, 3)
        P[o] = Po
    Rsp = minimal_measurements(F, n, X, P, groups, alpha)
    return P, groups, alpha, Rsp
