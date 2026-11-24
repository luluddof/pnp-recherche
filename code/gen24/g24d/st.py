# G24d: improved SAT encoding of the monomial ST model (sign-free relaxation + CEGAR for signs).
# Colours 0..q-1 = S, q..q+p-1 = T  (AT MOST q S-colours, p T-colours; unused allowed).
# Symmetry breaking (soundness argument in notes.md):
#  (B1) row 0 = a line (row or column, transpose allowed) maximising max(#S,#T) [balanced mode, p==q, palettes swapped so S is
#       the majority] or maximising #S [unbalanced mode]; s := #S in row 0. Every row/col then has #S<=s (and #T<=s if balanced).
#  (B2) columns: row 0 = S_0..S_{s-1} | T_0..T_{N-s-1}; inside each block the columns are sorted by #S (non-increasing).
#  (B3) rows 1..N-1 sorted by the colour index of column 0 (strictly increasing).
#  (B4) "new" colours (not in row 0) of each palette appear in first-occurrence order along col 0 (rows 1..), then rows 1.. .
#  (K)  class-size caps from the key lemma: |S-class| <= g*(p), |T-class| <= g*(q).
# Signs: the sign condition (product -1 on each intercalate) is a GF(2) linear system given the colouring; solved after each
# model; if inconsistent, a minimal-ish dependent set of intercalates is blocked (CEGAR). Sign-free UNSAT => UNSAT.
import sys, itertools, time, json, os
from pysat.solvers import Solver
from pysat.card import CardEnc, EncType

def gstar(N, t):
    g = 0
    while g + 1 <= N and (g + 1) * g // 2 <= t: g += 1
    return g

class Enc:
    def __init__(self): self.nv = 0; self.cl = []
    def new(self): self.nv += 1; return self.nv
    def add(self, c): self.cl.append(c)
    def amo(self, lits):
        if len(lits) <= 6:
            for a, b in itertools.combinations(lits, 2): self.add([-a, -b])
        else:
            e = CardEnc.atmost(lits, bound=1, top_id=self.nv, encoding=EncType.seqcounter)
            self.nv = max(self.nv, e.nv); self.cl += e.clauses
    def atmost(self, lits, k):
        if k >= len(lits): return
        e = CardEnc.atmost(lits, bound=k, top_id=self.nv, encoding=EncType.seqcounter)
        self.nv = max(self.nv, e.nv); self.cl += e.clauses
    def unary(self, lits):
        """full totalizer: returns o[0..n-1], o[t] <-> (sum >= t+1)."""
        if len(lits) == 1: return lits[:]
        h = len(lits) // 2
        A = self.unary(lits[:h]); B = self.unary(lits[h:])
        o = [self.new() for _ in range(len(lits))]
        for i in range(len(A) + 1):
            for j in range(len(B) + 1):
                if i + j >= 1:   # A>=i & B>=j -> o>=i+j
                    c = [o[i + j - 1]]
                    if i: c.append(-A[i - 1])
                    if j: c.append(-B[j - 1])
                    self.add(c)
                if i + j < len(o):  # A<i+1 & B<j+1 -> o < i+j+1
                    c = [-o[i + j]]
                    if i < len(A): c.append(A[i])
                    if j < len(B): c.append(B[j])
                    self.add(c)
        return o

def build(N, p, q, s, balanced=True, caps=True, colsort=True):
    E = Enc(); n = q + p
    isS = lambda k: k < q
    x = [[[E.new() for k in range(n)] for b in range(N)] for a in range(N)]
    P = [[E.new() for b in range(N)] for a in range(N)]
    # row 0
    r0 = [b if b < s else q + b - s for b in range(N)]
    assert s <= q and N - s <= p
    for b in range(N): E.add([x[0][b][r0[b]]])
    for a in range(N):
        for b in range(N):
            E.add(x[a][b][:]); E.amo(x[a][b])
            for k in range(n): E.add([-x[a][b][k], P[a][b] if isS(k) else -P[a][b]])
    for k in range(n):
        for a in range(N):
            E.amo([x[a][b][k] for b in range(N)]); E.amo([x[b][a][k] for b in range(N)])
    # S-colours of row 0 other than S_0 never appear in column 0 (implied; helps)
    for a in range(1, N):
        for k in range(1, s): E.add([-x[a][0][k]])
    # equality vars
    eqv = {}
    cells = [(a, b) for a in range(N) for b in range(N)]
    def eq(c1, c2):
        key = (min(c1, c2), max(c1, c2))
        if key in eqv: return eqv[key]
        v = E.new(); eqv[key] = v
        (a1, b1), (a2, b2) = key
        for k in range(n):
            E.add([-v, -x[a1][b1][k], x[a2][b2][k]])
            E.add([v, -x[a1][b1][k], -x[a2][b2][k]])
        return v
    inter = {}
    for r1, r2 in itertools.combinations(range(N), 2):
        for c1, c2 in itertools.combinations(range(N), 2):
            eD = eq((r1, c1), (r2, c2)); eA = eq((r1, c2), (r2, c1))
            E.add([-eD, eA]); E.add([eD, -eA])
            E.add([-eD, P[r1][c1], P[r1][c2]]); E.add([-eD, -P[r1][c1], -P[r1][c2]])
            inter[(r1, r2, c1, c2)] = eD
    # line counts
    rowc = [E.unary([P[a][b] for b in range(N)]) for a in range(N)]
    colc = [E.unary([P[a][b] for a in range(N)]) for b in range(N)]
    for o in rowc + colc:
        if s < N: E.add([-o[s]])                      # #S <= s
        if balanced and N - s - 1 >= 0: E.add([o[N - s - 1]])   # #T <= s  <=> #S >= N-s
    c0 = os.environ.get('C0')    # cube: column 0 has exactly C0 S-cells
    if c0 is not None:
        c0 = int(c0)
        if c0 >= 1: E.add([colc[0][c0 - 1]])
        if c0 < N: E.add([-colc[0][c0]])
    for item in filter(None, os.environ.get('CC', '').split(';')):   # cube: column b has exactly c S-cells
        b, c = map(int, item.split(':'))
        if c >= 1: E.add([colc[b][c - 1]])
        if c < N: E.add([-colc[b][c]])
    if colsort:
        for lo, hi in ((0, s), (s, N)):
            for b in range(lo, hi - 1):
                for t in range(N): E.add([-colc[b + 1][t], colc[b][t]])
    # column 0 rows 1..N-1 strictly increasing colour index
    ge = [[None] * (n + 1) for a in range(N)]
    for a in range(1, N):
        ge[a][n] = None
        prev = None
        for k in range(n - 1, -1, -1):
            v = E.new(); ge[a][k] = v
            E.add([-x[a][0][k], v]);
            if prev: E.add([-prev, v]); E.add([-v, x[a][0][k], prev])
            else: E.add([-v, x[a][0][k]])
            prev = v
    for a in range(1, N - 1):
        for k in range(n):
            # x[a+1][0][k] -> index(a,0) < k  i.e. not ge[a][k]
            E.add([-x[a + 1][0][k], -ge[a][k]])
    # first-occurrence order of new colours per palette
    scan = [(a, 0) for a in range(1, N)] + [(a, b) for a in range(1, N) for b in range(1, N)]
    for ks in (list(range(s, q)), list(range(q + N - s, n))):
        if len(ks) < 2: continue
        prev = None
        for i, (a, b) in enumerate(scan):
            cur = {}
            for k in ks:
                v = E.new(); cur[k] = v
                if prev is None: E.add([-v, x[a][b][k]]); E.add([v, -x[a][b][k]])
                else: E.add([-v, x[a][b][k], prev[k]]); E.add([v, -x[a][b][k]]); E.add([v, -prev[k]])
            for j in range(1, len(ks)):
                if prev is None: E.add([-x[a][b][ks[j]]])
                else: E.add([-x[a][b][ks[j]], prev[ks[j - 1]]])
            prev = cur
    if caps:
        cS, cT = gstar(N, p), gstar(N, q)
        for k in range(n):
            E.atmost([x[a][b][k] for a in range(N) for b in range(N)], cS if isS(k) else cT)
    lo = int(os.environ.get('CLASSMIN', '0'))   # SAT-HUNTING ONLY: every colour class has >= lo cells
    if lo:
        for k in range(n):
            e = CardEnc.atleast([x[a][b][k] for a in range(N) for b in range(N)], bound=lo, top_id=E.nv, encoding=EncType.seqcounter)
            E.nv = max(E.nv, e.nv); E.cl += e.clauses
    return E, x, inter

def gf2_signs(N, col):
    """solve s(a,b) in GF(2) with s11+s22+s12+s21 = 1 on every intercalate. returns (signs, None) or (None, conflict set)."""
    eqs = []; tags = []
    for r1, r2 in itertools.combinations(range(N), 2):
        for c1, c2 in itertools.combinations(range(N), 2):
            if col[r1][c1] == col[r2][c2]:
                m = (1 << (r1 * N + c1)) | (1 << (r2 * N + c2)) | (1 << (r1 * N + c2)) | (1 << (r2 * N + c1))
                eqs.append([m, 1, 1 << len(tags)]); tags.append((r1, r2, c1, c2))
    piv = {}
    for e in eqs:
        m, r, h = e
        while m:
            t = m.bit_length() - 1
            if t in piv:
                pm, pr, ph = piv[t]; m ^= pm; r ^= pr; h ^= ph
            else: piv[t] = (m, r, h); break
        if m == 0 and r == 1:
            return None, [tags[i] for i in range(len(tags)) if h >> i & 1]
    sol = 0
    for t in sorted(piv):
        m, r, h = piv[t]
        rest = m & ~(1 << t)
        val = r ^ (bin(rest & sol).count('1') & 1)
        if val: sol |= 1 << t
    sg = [[-1 if sol >> (a * N + b) & 1 else 1 for b in range(N)] for a in range(N)]
    return sg, None

def solve(N, p, q, s, balanced=True, solver='cadical195', tlimit=None, dimacs=None, verbose=True, **kw):
    t0 = time.time()
    E, x, inter = build(N, p, q, s, balanced, **kw)
    if dimacs:
        with open(dimacs, 'w') as f:
            f.write(f"p cnf {E.nv} {len(E.cl)}\n")
            for c in E.cl: f.write(' '.join(map(str, c)) + ' 0\n')
    if verbose: print(f"  enc N={N} p={p} q={q} s={s}: vars {E.nv} clauses {len(E.cl)} [{time.time()-t0:.1f}s]", flush=True)
    S = Solver(name=solver, bootstrap_with=E.cl)
    n = p + q; it = 0
    while True:
        it += 1
        r = S.solve()
        if not r: return 'UNSAT', None, it
        m = S.get_model(); ms = set(l for l in m if l > 0)
        col = [[next(k for k in range(n) if x[a][b][k] in ms) for b in range(N)] for a in range(N)]
        sg, conf = gf2_signs(N, col)
        if sg is not None: return 'SAT', (col, sg), it
        S.add_clause([-inter[t] for t in conf])
        if verbose: print(f"  sign conflict of size {len(conf)} blocked", flush=True)

def verify(N, p, q, col, sg):
    from core import predicate, ranks, identity_tensor_ok
    pal = {k: ('S' if k < q else 'T') for k in range(p + q)}
    used = set(c for r in col for c in r)
    nS = sum(1 for k in used if k < q); nT = len(used) - nS
    pr = predicate(N, col, pal, sg)
    idt = identity_tensor_ok(N, col, pal, sg)
    ok, rr, rc = ranks(N, col, pal, sg)
    # independent check with the g22c evaluator: C = sum_S K (x) K^T, p = rank(X - C)
    import numpy as np, ev
    terms = []
    for k in used:
        if k < q:
            K = np.zeros((N, N), dtype=np.int64)
            for a in range(N):
                for b in range(N):
                    if col[a][b] == k: K[a, b] = sg[a][b]
            terms.append((K, K.T))
    pe, meth = ev.p_of(N, terms)
    return dict(nS=nS, nT=nT, predicate=pr, identity=idt, RplusC=ok, rankR=rr, rankC=rc, ev_rankR=pe, ev_meth=meth,
                ok=bool(pr and idt and ok and rr <= p and rc <= q and pe <= p))

if __name__ == "__main__":
    N, p, q = map(int, sys.argv[1:4])
    ss = [int(t) for t in sys.argv[4].split(',')] if len(sys.argv) > 4 and sys.argv[4] != 'all' else None
    bal = (p == q) and not os.environ.get('NOBAL')
    kw = dict(colsort=not os.environ.get('NOCOLSORT'), caps=not os.environ.get('NOCAPS'))
    solver = os.environ.get('SOLVER', 'cadical195')
    if ss is None:
        ss = [s for s in range(N + 1) if s <= q and N - s <= p and (not bal or 2 * s >= N)]
    tag = os.environ.get('TAG', '')
    for s in ss:
        t = time.time()
        st, sol, it = solve(N, p, q, s, balanced=bal, solver=solver, dimacs=os.environ.get('DIMACS', '').replace('@', str(s)) or None, **kw)
        msg = f"N={N} (p,q)=({p},{q}) s={s}: {st} it={it} [{time.time()-t:.1f}s]"
        if sol:
            col, sg = sol
            v = verify(N, p, q, col, sg); msg += ' ' + json.dumps(v)
            fn = f"sol_{N}_{p}_{q}_s{s}.json"
            json.dump(dict(N=N, p=p, q=q, s=s, col=col, sg=sg, pal='k<q:S else T', verify=v), open(fn, 'w'))
            print(msg, flush=True); break
        print(msg, flush=True)
