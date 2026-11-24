"""G14c library: interface-rank (IR) model for clause-product schemes. Exact arithmetic (Fractions)."""
import sys, random, itertools, math
sys.dont_write_bytecode = True
sys.path.insert(0, '/home/user/pnp-recherche/code/gen6/g6a')   # read-only import
from fractions import Fraction as Fr
from scheme import rank_Q, Pval, occurrences, random_3regular_formula
import designs as D

def form_matrix(F, n, P, t):
    """A_t: dict occurrence -> row (list of n Fractions) = P_o(t)."""
    A = {}
    for o in occurrences(F):
        val = Pval(P[o], t)
        A[o] = [val.get(u, Fr(0)) for u in range(n)]
    return A

def trivial_P(F, n):
    return {o: {o[1]: [Fr(1)]} for o in occurrences(F)}

class Cuts:
    """rank tables r[Emask][Smask] = rank A[rows of clauses in Emask, columns Smask]; kappa on leaf masks.
    Leaves: clause e -> bit e (e<m), variable v -> bit m+v."""
    def __init__(self, F, n, A):
        self.F, self.n, self.m, self.A = F, n, len(F), A
        m = self.m
        self.rows_of = [[A[(e, v)] for (v, _) in F[e]] for e in range(m)]
        self.r = {}
        self.fullE, self.fullS = (1 << m) - 1, (1 << n) - 1
    def rk(self, Em, Sm):
        key = (Em, Sm)
        if key not in self.r:
            cols = [u for u in range(self.n) if Sm >> u & 1]
            rows = [[row[u] for u in cols] for e in range(self.m) if Em >> e & 1 for row in self.rows_of[e]]
            self.r[key] = rank_Q(rows) if cols and rows else 0
        return self.r[key]
    def kappa_ES(self, Em, Sm):
        return self.rk(Em, self.fullS ^ Sm) + self.rk(self.fullE ^ Em, Sm)
    def kappa(self, mask):
        return self.kappa_ES(mask & self.fullE, mask >> self.m)

def branchwidth(C):
    """exact min over cubic trees on the m+n leaves of max edge kappa (subset DP, 3^L)."""
    L = C.m + C.n; full = (1 << L) - 1
    kap = [0] * (1 << L)
    for S in range(1, full):
        kap[S] = C.kappa(S)
    f = [0] * (1 << L)
    order = sorted(range(1, full + 1), key=lambda S: bin(S).count('1'))
    for S in order:
        if S & (S - 1) == 0:
            f[S] = kap[S]; continue
        low = S & -S; rest = S ^ low; best = 10 ** 9
        sub = rest
        while True:            # S1 = low | (subset of rest), S2 = S ^ S1 nonempty
            S2 = rest ^ sub    # sub = part of rest joining low ; S2 = complement in rest
            if S2:
                v = max(f[low | sub], f[S2])
                if v < best: best = v
            if sub == 0: break
            sub = (sub - 1) & rest
        f[S] = max(kap[S], best) if S != full else best
    return f[full]

def treewidth_exact(N, adj):
    """exact treewidth (width) by subset DP over elimination prefixes; adj: dict v->set. N<=16."""
    if N == 0: return -1
    INF = 10 ** 9
    TW = {0: -1}
    def q(Smask, v):
        seen = 1 << v; stack = [v]; cnt = 0
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if seen >> y & 1: continue
                seen |= 1 << y
                if Smask >> y & 1: stack.append(y)
                else: cnt += 1
        return cnt
    full = (1 << N) - 1
    tw = [INF] * (1 << N); tw[0] = -1
    for S in range(1, full + 1):
        best = INF
        for v in range(N):
            if S >> v & 1:
                T = S ^ (1 << v)
                c = max(tw[T], q(T, v))
                if c < best: best = c
        tw[S] = best
    return tw[full]

def primal_adj(F, n):
    adj = {v: set() for v in range(n)}
    for c in F:
        vs = [v for v, _ in c]
        for a in vs:
            for b in vs:
                if a != b: adj[a].add(b)
    return adj

def coords(basis, vec):
    """vec = sum lam_j basis[j] (basis independent, vec in span): return lam."""
    k = len(basis)
    if k == 0:
        assert not any(vec); return []
    ncol = len(vec)
    M = [[basis[j][i] for j in range(k)] + [vec[i]] for i in range(ncol)]
    piv = []; r = 0
    for c in range(k):
        p = next((i for i in range(r, ncol) if M[i][c] != 0), None)
        assert p is not None
        M[r], M[p] = M[p], M[r]
        pv = M[r][c]; M[r] = [x / pv for x in M[r]]
        for i in range(ncol):
            if i != r and M[i][c] != 0:
                fct = M[i][c]; M[i] = [a - fct * b for a, b in zip(M[i], M[r])]
        piv.append(r); r += 1
    assert all(M[i][k] == 0 for i in range(r, ncol)), "vec not in span"
    return [M[j][k] for j in range(k)]

def row_graph(F, n, A, groups, alpha, Rsp):
    """Row graph in the basis Rsp[h] (effective spaces R'_h): vertices (h,j); for each group g a clique on
    its own rows and the basis rows of other groups in the support of its foreign restrictions."""
    V = [(h, j) for h in range(len(groups)) for j in range(len(Rsp[h]))]
    ix = {v: i for i, v in enumerate(V)}
    scope = [set(ix[(g, j)] for j in range(len(Rsp[g]))) for g in range(len(groups))]
    reads = [set() for _ in groups]           # groups read by g  (L_{g->h} != 0)
    for (e, v), row in A.items():
        g = alpha[e]
        for h, S in enumerate(groups):
            if h == g: continue
            vec = [row[u] for u in S]
            if any(vec):
                reads[g].add(h)
                lam = coords(Rsp[h], vec)
                for j, l in enumerate(lam):
                    if l != 0: scope[g].add(ix[(h, j)])
    adj = {i: set() for i in range(len(V))}
    for sc in scope:
        for a in sc:
            for b in sc:
                if a != b: adj[a].add(b)
    return V, adj, scope, reads

