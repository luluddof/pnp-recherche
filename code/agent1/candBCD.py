"""Candidates B (LP: cluster/Sherali-Adams-style local-distribution LP + LP rounding
+ poly-budget WalkSAT), C (BP-guided decimation + unit propagation),
D (spectral refutation, Feige-Ofek / Coja-Oghlan-Goerdt-Lanka style)."""
import itertools, math, random
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog
from scipy.sparse.linalg import eigsh


# ---------------- fast unit propagation ----------------
def make_occ(n, cls):
    occ = [[] for _ in range(2 * n + 2)]
    for ci, c in enumerate(cls):
        for l in c:
            occ[l].append(ci)
    return occ


def up(n, cls, occ, assign, queue):
    """assign: list (None/True/False) of len n+1, modified in place. queue: literals newly set true."""
    q = list(queue)
    while q:
        l = q.pop()
        for ci in occ[-l]:          # clauses where -l occurs: may become unit/empty
            unk = None; nunk = 0; sat = False
            for x in cls[ci]:
                v = assign[abs(x)]
                if v is None:
                    nunk += 1; unk = x
                elif v == (x > 0):
                    sat = True; break
            if sat:
                continue
            if nunk == 0:
                return False
            if nunk == 1:
                assign[abs(unk)] = unk > 0
                q.append(unk)
    return True


# ---------------- B: LPs ----------------
def cluster_lp_feasible(n, cls):
    """Local-distribution LP: a distribution over satisfying assignments of each clause
    and of each union of two clauses that share a variable; marginals consistent on
    every clause scope. Infeasible => UNSAT (sound). Size poly(m) (clusters <= m^2, each <= 2^6)."""
    scopes = [tuple(sorted(set(abs(l) for l in c))) for c in cls]
    var_of = {}
    nv = 0
    def newvar():
        nonlocal nv; nv += 1; return nv - 1
    rows, cols, vals, rhs = [], [], [], []
    r = 0
    def satisfies(c, scope, a):
        pos = {v: i for i, v in enumerate(scope)}
        return any(a[pos[abs(l)]] == (l > 0) for l in c)
    # clause-level distributions q_C
    qvars = []
    for ci, c in enumerate(cls):
        sc = scopes[ci]
        d = {}
        for a in itertools.product([False, True], repeat=len(sc)):
            if satisfies(c, sc, a):
                d[a] = newvar()
        qvars.append(d)
        for vid in d.values():
            rows.append(r); cols.append(vid); vals.append(1.)
        rhs.append(1.); r += 1
    occ = {}
    for ci, sc in enumerate(scopes):
        for v in sc:
            occ.setdefault(v, []).append(ci)
    # singleton consistency x_v across clauses: marginal of v in q_C equal across C
    for v, cl in occ.items():
        for c1, c2 in zip(cl, cl[1:]):
            for C, sgn in ((c1, 1.), (c2, -1.)):
                sc = scopes[C]; p = sc.index(v)
                for a, vid in qvars[C].items():
                    if a[p]:
                        rows.append(r); cols.append(vid); vals.append(sgn)
            rhs.append(0.); r += 1
    pairs = set()
    for v, cl in occ.items():
        for c1, c2 in itertools.combinations(cl, 2):
            pairs.add((c1, c2))
    for c1, c2 in pairs:
        K = tuple(sorted(set(scopes[c1]) | set(scopes[c2])))
        if len(K) > 6:
            continue
        pv = {}
        for a in itertools.product([False, True], repeat=len(K)):
            if satisfies(cls[c1], K, a) and satisfies(cls[c2], K, a):
                pv[a] = newvar()
        for C in (c1, c2):
            sc = scopes[C]; idx = [K.index(v) for v in sc]
            for b, qid in qvars[C].items():
                for a, vid in pv.items():
                    if tuple(a[i] for i in idx) == b:
                        rows.append(r); cols.append(vid); vals.append(1.)
                rows.append(r); cols.append(qid); vals.append(-1.)
                rhs.append(0.); r += 1
            # assignments of sc that falsify C have no q var but must also get 0 mass: implied by pv support
    A = sp.csr_matrix((vals, (rows, cols)), shape=(r, nv))
    res = linprog(np.zeros(nv), A_eq=A, b_eq=np.array(rhs), bounds=(0, 1), method='highs')
    return res.status == 0   # 0 feasible, 2 infeasible


def basic_lp(n, cls, rng):
    """Standard LP: sum of literal values >= 1 per clause, random objective."""
    A = []; b = []
    rows, cols, vals = [], [], []
    for i, c in enumerate(cls):
        const = 0
        for l in c:
            if l > 0:
                rows.append(i); cols.append(l - 1); vals.append(-1.)
            else:
                rows.append(i); cols.append(-l - 1); vals.append(1.); const += 1
        b.append(-1. + const)
    A = sp.csr_matrix((vals, (rows, cols)), shape=(len(cls), n))
    obj = np.array([rng.uniform(-1, 1) for _ in range(n)])
    res = linprog(obj, A_ub=A, b_ub=np.array(b), bounds=(0, 1), method='highs')
    return res.x if res.status == 0 else None


def walksat(n, cls, init, max_flips, rng, p=0.5):
    a = list(init)
    occ = make_occ(n, cls)
    nsat = [sum(1 for l in c if a[abs(l)] == (l > 0)) for c in cls]
    unsat = set(i for i, k in enumerate(nsat) if k == 0)
    for _ in range(max_flips):
        if not unsat:
            return a
        ci = rng.choice(tuple(unsat)) if len(unsat) < 64 else rng.sample(sorted(unsat), 1)[0]
        c = cls[ci]
        def breaks(v):
            lit = v if a[v] else -v   # literal currently true
            return sum(1 for cj in occ[lit] if nsat[cj] == 1)
        if rng.random() < p:
            v = abs(rng.choice(c))
        else:
            v = min((abs(l) for l in c), key=breaks)
        lit = v if a[v] else -v
        for cj in occ[lit]:
            nsat[cj] -= 1
            if nsat[cj] == 0: unsat.add(cj)
        a[v] = not a[v]
        for cj in occ[-lit]:
            nsat[cj] += 1
            unsat.discard(cj)
    return None


def candB(n, cls, rng, flips_factor=20):
    """SAT side: LP solution, randomized rounding, then WalkSAT with flips_factor*n^2 flips."""
    x = basic_lp(n, cls, rng)
    init = [None] + [rng.random() < x[i] for i in range(n)]
    a = walksat(n, cls, init, flips_factor * n * n, rng)
    return a


# ---------------- C: BP-guided decimation ----------------
def bp_decimate(n, cls, rng, iters=60, frac=0.05, damp=0.5):
    """Polynomial: <= n rounds, each <= iters*O(m) BP work + O(nm) unit propagation."""
    assign = [None] * (n + 1)
    occ = make_occ(n, cls)
    while True:
        live = []; lits = []
        for c in cls:
            if any(assign[abs(l)] == (l > 0) for l in c if assign[abs(l)] is not None):
                continue
            rest = [l for l in c if assign[abs(l)] is None]
            if not rest:
                return None
            live.append(rest)
        free = [v for v in range(1, n + 1) if assign[v] is None]
        if not free:
            return assign
        if not live:
            for v in free: assign[v] = True
            return assign
        E_c = np.array([ci for ci, c in enumerate(live) for l in c])
        E_v = np.array([abs(l) for c in live for l in c])
        E_s = np.array([l > 0 for c in live for l in c])
        E = len(E_c)
        eta = np.full(E, 0.5) * rng.random() if False else np.array([rng.random() for _ in range(E)])
        for _ in range(iters):
            lg = np.log(np.clip(1 - eta, 1e-300, 1))
            Spos = np.bincount(E_v[E_s], weights=lg[E_s], minlength=n + 1)
            Sneg = np.bincount(E_v[~E_s], weights=lg[~E_s], minlength=n + 1)
            same = np.where(E_s, Spos[E_v], Sneg[E_v]) - lg
            opp = np.where(E_s, Sneg[E_v], Spos[E_v])
            # pi = Pu/(Pu+Ps), Pu=exp(same), Ps=exp(opp)
            pi = 1 / (1 + np.exp(np.clip(opp - same, -700, 700)))
            lpi = np.log(np.clip(pi, 1e-300, 1))
            tot = np.bincount(E_c, weights=lpi, minlength=len(live))
            new = np.exp(tot[E_c] - lpi)
            eta = damp * eta + (1 - damp) * new
        lg = np.log(np.clip(1 - eta, 1e-300, 1))
        Spos = np.bincount(E_v[E_s], weights=lg[E_s], minlength=n + 1)
        Sneg = np.bincount(E_v[~E_s], weights=lg[~E_s], minlength=n + 1)
        # P(x=T) ∝ prod over clauses where v negative (1-eta) = exp(Sneg)
        pT = 1 / (1 + np.exp(np.clip(Spos - Sneg, -700, 700)))
        order = sorted(free, key=lambda v: -abs(pT[v] - 0.5))
        k = max(1, int(frac * len(free)))
        for v in order[:k]:
            if assign[v] is not None:
                continue
            val = bool(pT[v] >= 0.5)
            assign[v] = val
            if not up(n, cls, occ, assign, [v if val else -v]):
                return None


# ---------------- D: spectral refutation ----------------
def spectral_refute(n, cls):
    """Sound certificate: for every x in {±1}^n, #violated = (1/8) sum_C (1 - e1 + e2 - e3) and
    sum_C(e1 - e2 + e3) <= L + n*lmax(A) + sqrt(n (tr G + n^2 lmax(G_off))). If bound < m => UNSAT.
    Polytime: one eigenvalue of an n x n and an n^2 x n^2 sparse matrix."""
    m = len(cls)
    c = np.zeros(n + 1)
    A = np.zeros((n + 1, n + 1))
    trip = []
    for cl in cls:
        s = [1 if l > 0 else -1 for l in cl]; v = [abs(l) for l in cl]
        for i in range(3):
            c[v[i]] += s[i]
        for i, j in ((0, 1), (0, 2), (1, 2)):
            A[v[i], v[j]] -= s[i] * s[j] / 2; A[v[j], v[i]] -= s[i] * s[j] / 2
        trip.append((v, s[0] * s[1] * s[2]))
    L = np.abs(c).sum()
    Qb = n * max(0., np.linalg.eigvalsh(A[1:, 1:])[-1])
    # B_i[j,k] entries: s/6 for all 6 orderings
    rows, cols, vals = [], [], []  # W: n x n^2, row i = vec(B_i)
    for v, s in trip:
        for i, j, k in itertools.permutations(range(3)):
            rows.append(v[i] - 1); cols.append((v[j] - 1) * n + (v[k] - 1)); vals.append(s / 6)
    W = sp.csr_matrix((vals, (rows, cols)), shape=(n, n * n))
    G = (W.T @ W).tocsr()
    trG = G.diagonal().sum()
    Goff = G - sp.diags(G.diagonal())
    try:
        lg = eigsh(Goff.astype(float), k=1, which='LA', return_eigenvectors=False)[0]
    except Exception:
        lg = np.linalg.eigvalsh(Goff.toarray())[-1]
    Tb = math.sqrt(n * (trG + n * n * max(0., lg)))
    bound = L + Qb + Tb
    return bound < m, bound / m
