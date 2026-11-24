# One-copy objects for the flattening bound:  B_g^{(ha,hb)} (pairs) and O_g^{(e)} (labels), exact rationals.
#   F^{(e)}[ha,hb] = sum_g B_g^{(ha,hb)} (x) O_g^{(e)}     rows (a, m'), cols (b, m)      ("cross flattening 1")
#   F2             = sum_g B_g (x) O_g^T                    rows (a, m),  cols (b, m')     ("cross flattening 2")
# For ANY K=0 construction: target delta_{xv} delta_{uy} flattened (x,u)|(v,y) is I_{N^2}, so N^2 <= rank F (and F2).
import sys
sys.path.insert(0, '/tmp/claude-0/-home-user-pnp-recherche/2d12e629-e6d2-542b-bfc0-8b766f57383f/scratchpad/gen24/g24f')
from trank import *
from fractions import Fraction as Fr
from functools import lru_cache

class Copy:
    def __init__(self, G, xi, coset=None):
        self.G = G; self.xi = tuple(Fr(x) for x in xi); self.r = len(G); self.coset = coset
        self.S = [tuple(int(c) for c in g) for g in S_xi(G, self.xi)]
        self._basis = {}
    def basis(self, h, module=False):
        key = (h, module)
        if key not in self._basis:
            cos = self.coset if module else None
            out = []
            for p in lattice_points(self.G, 2*(h + 2*ip(self.G, self.xi, self.xi) + 4), cos):
                k = h - wt(self.G, self.xi, p)
                if k >= 0 and k.denominator == 1:
                    for m in fock(self.r, int(k)): out.append((p, m))
            self._basis[key] = out
        return self._basis[key]
    def B(self, g, ha, hb):
        A = self.basis(ha); Bs = self.basis(hb); r = self.r
        M = [[Fr(0)]*len(Bs) for _ in A]
        for i, (beta, am) in enumerate(A):
            for j, (delta, bm) in enumerate(Bs):
                if all(beta[c]+delta[c] == g[c] for c in range(r)):
                    M[i][j] = Bval(self.G, beta, am, delta, bm)
        return M
    def O(self, g, e):
        Mb = self.basis(e, module=True); idx = {b: i for i, b in enumerate(Mb)}
        n = len(Mb); M = [[Fr(0)]*n for _ in range(n)]   # M[out][in]
        for j, (p, mono) in enumerate(Mb):
            q = tuple(p[c]+g[c] for c in range(self.r)); k2 = e - wt(self.G, self.xi, q)
            if k2 < 0 or k2.denominator != 1: continue
            for m2, c2 in O_apply(self.G, g, p, mono, int(k2)).items():
                M[idx[(q, m2)]][j] += c2
        return M

def kron(A, B):
    return [[a*b for a in ra for b in rb] for ra in A for rb in B]
def madd(A, B):
    return [[x+y for x, y in zip(ra, rb)] for ra, rb in zip(A, B)]
def transpose(A):
    return [list(r) for r in zip(*A)] if A else A
def rankQ(M):
    if not M or not M[0]: return 0
    # clear denominators row-wise, exact integer rank
    import math
    rows = []
    for r in M:
        d = 1
        for x in r: d = d*x.denominator//math.gcd(d, x.denominator)
        rows.append([int(x*d) for x in r])
    return flint.fmpz_mat(rows).rank()

def flat_rank(C, e, ha, hb, which=1, return_blocks=False):
    """exact rank of F = sum_g B_g (x) O_g (which=1, rows (a,m') cols (b,m)) or sum_g B_g (x) O_g^T (which=2,
    rows (a,m) cols (b,m')), computed per conserved-charge block."""
    import math
    A = C.basis(ha); Bb = C.basis(hb); Mb = C.basis(e, True)
    Bs = {g: C.B(g, ha, hb) for g in C.S}; Os = {g: C.O(g, e) for g in C.S}
    r = C.r
    # rows (i,k), cols (j,l): entry sum_g Bs[g][i][j] * (Os[g][k][l] if which==1 else Os[g][l][k])
    # conserved: which=1: charge(m') - charge(a) = charge(m) + charge(b); which=2: charge(m)+charge(a)... use generic grouping by
    # connected components instead (robust): build sparse entries then group.
    from collections import defaultdict
    ent = defaultdict(Fr)
    nzB = {g: [(i, j, Bs[g][i][j]) for i in range(len(A)) for j in range(len(Bb)) if Bs[g][i][j]] for g in C.S}
    nzO = {g: [(k, l, Os[g][k][l]) for k in range(len(Mb)) for l in range(len(Mb)) if Os[g][k][l]] for g in C.S}
    for g in C.S:
        for (i, j, x) in nzB[g]:
            for (k, l, y) in nzO[g]:
                if which == 1: ent[((i, k), (j, l))] += x*y
                else: ent[((i, l), (j, k))] += x*y
    # union-find on rows/cols
    par = {}
    def f(u):
        while par.setdefault(u, u) != u:
            par[u] = par[par[u]]; u = par[u]
        return u
    for (ro, co), v in ent.items():
        if v: par[f(('r', ro))] = f(('c', co))
    comps = defaultdict(lambda: (set(), set()))
    for (ro, co), v in ent.items():
        if v:
            c = f(('r', ro)); comps[c][0].add(ro); comps[c][1].add(co)
    tot = 0; blocks = []
    for c, (rs, cs) in comps.items():
        rs = sorted(rs); cs = sorted(cs); ri = {x: n for n, x in enumerate(rs)}; ci = {x: n for n, x in enumerate(cs)}
        M = [[Fr(0)]*len(cs) for _ in rs]
        for (ro, co), v in ent.items():
            if v and ro in ri and co in ci: M[ri[ro]][ci[co]] = v
        rk = rankQ(M); tot += rk; blocks.append((len(rs), len(cs), rk))
    return (tot, blocks) if return_blocks else tot
