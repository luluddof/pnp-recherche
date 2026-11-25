# G26-P4: exact search in the monomial-C model MC(N,p,q).
# K[(x,y),(u,v)] = sum_l E_l[x,u] C_l[y,v], E_l, C_l signed partial permutations; need rank_Q(Swap - K) <= p.
# Stage 1 (sound relaxation): GF(2) SAT on supports: rank_GF2(Swap + K mod 2) <= p via M = U V, V in RREF.
#   GF(2)-UNSAT (for every normal-form branch r) => MC-UNSAT, because rank_Q(M) >= rank_GF2(M mod 2).
# Symmetry breaking (valid in MC): group (A,B) perms: E_l -> A E_l B^T, C_l -> B C_l A^T; term order; flip E<->C.
#   r = max rank over all 2q matrices; WLOG it is E_1 and E_1 = diag(I_r,0); all |E_l|,|C_l| <= r;
#   C_1's lower-right (N-r)x(N-r) block in form diag(I_s,0) (independent perms of the complement);
#   |E_2| >= |E_3| >= ... (term order of l >= 2) [optional, ORD=1].
# Stage 2: for each GF(2) support, signs by GF(3) SAT (rank_GF3 <= p), then exact rank over Q (flint); CEGAR.
import sys, itertools, time, json, os
import numpy as np
from pysat.card import CardEnc, EncType, ITotalizer
from pysat.formula import IDPool

class Enc:
    def __init__(self): self.nv = 0; self.cl = []; self.xors = []
    def new(self): self.nv += 1; return self.nv
    def add(self, c): self.cl.append(c)

def rref(E, A, k, m):
    z = [[None] * (m + 1) for i in range(k)]
    PV = [[None] * m for i in range(k)]
    T = E.new(); E.add([T])
    for i in range(k):
        z[i][0] = T
        for j in range(m):
            v = E.new(); z[i][j + 1] = v
            E.add([-v, z[i][j]]); E.add([-v, -A[i][j]]); E.add([v, -z[i][j], A[i][j]])
    for i in range(k):
        for j in range(m):
            pv = E.new(); PV[i][j] = pv
            E.add([-pv, z[i][j]]); E.add([-pv, A[i][j]]); E.add([pv, -z[i][j], -A[i][j]])
            if i + 1 < k: E.add([-pv, z[i + 1][j + 1]])
            for i2 in range(k):
                if i2 != i: E.add([-pv, -A[i2][j]])
        if i + 1 < k: E.add([-z[i][m], z[i + 1][m]])
    return z, PV

def amo(E, lits):
    for a, b in itertools.combinations(lits, 2): E.add([-a, -b])

def unary(E, xs, K):
    """full-equivalence unary counter: returns o[k] <=> (sum xs >= k+1), k < K."""
    T = E.new(); E.add([T]); F = -T
    prev = [T] + [F] * K          # prev[k] <=> at least k among processed
    for x in xs:
        cur = [T]
        for k in range(1, K + 1):
            a, b = prev[k], prev[k - 1]
            z = E.new()           # z <=> a | (b & x)
            E.add([-a, z]); E.add([-b, -x, z]); E.add([-z, a, b]); E.add([-z, a, x])
            cur.append(z)
        prev = cur
    return prev[1:]

def card_le(E, lits, k):
    if k >= len(lits): return
    c = CardEnc.atmost(lits=lits, bound=k, top_id=E.nv, encoding=EncType.seqcounter)
    E.nv = max(E.nv, c.nv)
    for cl in c.clauses: E.add(cl)

def build(N, p, q, r, sym=True, order=False, extra=None, c22=True):
    E = Enc(); M = N * N
    e = [[[E.new() for u in range(N)] for x in range(N)] for l in range(q)]
    c = [[[E.new() for v in range(N)] for y in range(N)] for l in range(q)]
    for l in range(q):
        for X in (e[l], c[l]):
            for i in range(N):
                amo(E, [X[i][j] for j in range(N)]); amo(E, [X[j][i] for j in range(N)])
    U = [[E.new() for k in range(p)] for a in range(M)]
    V = [[E.new() for b in range(M)] for k in range(p)]
    def AND(a, b):
        z = E.new(); E.add([-z, a]); E.add([-z, b]); E.add([z, -a, -b]); return z
    LINK = os.environ.get('LINK', '1') == '1'
    Mv = {}
    NOUV = os.environ.get('NOUV') == '1'
    if NOUV: sym = False
    for x, y, u, v in (itertools.product(range(N), repeat=4) if not NOUV else []):
        terms = [AND(U[x * N + y][k], V[k][u * N + v]) for k in range(p)]
        kt = [AND(e[l][x][u], c[l][y][v]) for l in range(q)]
        if LINK:
            mv = E.new(); Mv[(x * N + y, u * N + v)] = mv
            E.xors.append((kt + [mv], int(x == v and y == u)))   # mv = Swap + K  (mod 2)
            E.xors.append((terms + [mv], 0))                     # mv = (U V)
        else:
            E.xors.append((terms + kt, int(x == v and y == u)))
    if r is not None:
        for x in range(N):
            for u in range(N): E.add([e[0][x][u]] if (x == u and x < r) else [-e[0][x][u]])
        for l in range(q):
            card_le(E, [e[l][i][j] for i in range(N) for j in range(N)], r)
            card_le(E, [c[l][i][j] for i in range(N) for j in range(N)], r)
        # C_1 lower-right block normal form diag(I_s,0)
        for i in (range(r, N) if c22 else []):
            for j in range(r, N):
                if i != j: E.add([-c[0][i][j]])
            if i + 1 < N: E.add([-c[0][i + 1][i + 1], c[0][i][i]])
        if order and q >= 3:
            outs = []
            for l in range(1, q):
                outs.append(unary(E, [e[l][i][j] for i in range(N) for j in range(N)], N))
            for a in range(len(outs) - 1):
                for k in range(N):
                    E.add([-outs[a + 1][k], outs[a][k]])
    if r is not None and os.environ.get('RSUM', '1') == '1':
        # implied constraint (rank subadditivity over any field): sum_l |E_l| |C_l| >= rank(K) >= N^2 - p
        prods = []
        for l in range(q):
            a = unary(E, [e[l][i][j] for i in range(N) for j in range(N)], r)
            b = unary(E, [c[l][i][j] for i in range(N) for j in range(N)], r)
            prods += [AND(x1, x2) for x1 in a for x2 in b]
        need = N * N - p
        if need > len(prods): E.add([])   # empty clause: branch impossible
        elif need > 0:
            cc = CardEnc.atleast(lits=prods, bound=need, top_id=E.nv, encoding=EncType.seqcounter)
            E.nv = max(E.nv, cc.nv)
            for x in cc.clauses: E.add(x)
    if os.environ.get('CORE', '1') == '1':
        # implied (first peeling round): Fix(Theta) embeds in span(dom ∩ im) => |dom ∩ im| >= N^2 - p, where
        # cell (a,b) in dom P_l iff col b of E_l and col a of C_l nonempty; in im P_l iff row a of E_l and row b of C_l nonempty
        def OR(lits):
            z = E.new(); E.add([-z] + lits)
            for t in lits: E.add([z, -t])
            return z
        rE = [[OR([e[l][i][j] for j in range(N)]) for i in range(N)] for l in range(q)]
        cE = [[OR([e[l][j][i] for j in range(N)]) for i in range(N)] for l in range(q)]
        rC = [[OR([c[l][i][j] for j in range(N)]) for i in range(N)] for l in range(q)]
        cC = [[OR([c[l][j][i] for j in range(N)]) for i in range(N)] for l in range(q)]
        both = []
        for a in range(N):
            for b in range(N):
                dd = OR([AND(cE[l][b], cC[l][a]) for l in range(q)])
                ii = OR([AND(rE[l][a], rC[l][b]) for l in range(q)])
                both.append(AND(dd, ii))
        need = N * N - p
        if need > 0:
            cc = CardEnc.atleast(lits=both, bound=need, top_id=E.nv, encoding=EncType.seqcounter)
            E.nv = max(E.nv, cc.nv)
            for x in cc.clauses: E.add(x)
    if os.environ.get('STAIR', '2') == '1':
        # [P] staircase lemma (MC form of g24g Lemma C, no cap): for any column order beta, rank >= sum_x t_beta(x),
        # t_beta(x) = number of leading columns (beta order) on which row x of the union support Z vanishes.
        # max over beta = DP over column subsets: h(S) = g(S) + max_{u in S} h(S - u), g(S) = #{x : Z[x,S] = 0}; need h <= p.
        def ORl(lits):
            z = E.new(); E.add([-z] + lits)
            for t in lits: E.add([z, -t])
            return z
        ZE = [[ORl([e[l][x][u] for l in range(q)]) for u in range(N)] for x in range(N)]
        ZC = [[ORl([c[l][x][u] for l in range(q)]) for u in range(N)] for x in range(N)]
        T = E.new(); E.add([T])
        for Z in (ZE, [list(r_) for r_ in zip(*ZE)], ZC, [list(r_) for r_ in zip(*ZC)]):
            w = {0: [T] * N}; H = {0: [T] + [None] * p}
            for S in range(1, 1 << N):
                hi = S.bit_length() - 1; Sp = S ^ (1 << hi)
                ws = []
                for x in range(N):
                    z = E.new(); E.add([-w[Sp][x], Z[x][hi], z]); ws.append(z)    # row x zero on S => z
                w[S] = ws
                G = unary(E, ws, min(N, p + 1))      # G[i-1] <=> g(S) >= i
                Gs = [T] + G
                Hs = [T] + [E.new() for j in range(p)]   # Hs[j] : h(S) >= j  (implied direction only)
                for u in range(N):
                    if not (S >> u) & 1: continue
                    Hp = H[S ^ (1 << u)]
                    for j in range(p + 1):
                        if Hp[j] is None: continue
                        for i in range(len(Gs)):
                            if i == 0 and j == 0: continue
                            if j + i >= p + 1: E.add([-Hp[j], -Gs[i]])
                            else: E.add([-Hp[j], -Gs[i], Hs[j + i]])
                H[S] = Hs
    if os.environ.get('STAIR', '2') == '2':
        # [P] combined bound (stair2.py validated): rank >= (N^2 - |D'|) + max_beta #{(x,y) in D' : pos(y) < t_beta(x)},
        # D' = dom ∩ im of Theta (support level), for the 4 variants (E,C), (E^T,C^T), (C,E), (C^T,E^T).
        def ORl(lits):
            z = E.new(); E.add([-z] + lits)
            for t in lits: E.add([z, -t])
            return z
        T = E.new(); E.add([T])
        tr = lambda L: [[[L[l][j][i] for j in range(N)] for i in range(N)] for l in range(q)]
        for EL, CL in ((e, c), (tr(e), tr(c)), (c, e), (tr(c), tr(e))):
            rowE = [[ORl(EL[l][i]) for i in range(N)] for l in range(q)]
            colE = [[ORl([EL[l][j][i] for j in range(N)]) for i in range(N)] for l in range(q)]
            rowC = [[ORl(CL[l][i]) for i in range(N)] for l in range(q)]
            colC = [[ORl([CL[l][j][i] for j in range(N)]) for i in range(N)] for l in range(q)]
            both = [[None] * N for _ in range(N)]
            for a in range(N):
                for b in range(N):
                    dd = ORl([AND(colE[l][b], colC[l][a]) for l in range(q)])
                    ii = ORl([AND(rowE[l][a], rowC[l][b]) for l in range(q)])
                    both[a][b] = AND(dd, ii)
            Z = [[ORl([EL[l][x][u] for l in range(q)]) for u in range(N)] for x in range(N)]
            Pc = unary(E, [-both[a][b] for a in range(N) for b in range(N)], p + 1)
            E.add([-Pc[p]])
            H = {0: [T] + Pc[:p]}
            w = {0: [T] * N}
            for S in range(1, 1 << N):
                hi = S.bit_length() - 1; Sp = S ^ (1 << hi)
                ws = []
                for x in range(N):
                    z = E.new(); E.add([-w[Sp][x], Z[x][hi], z]); ws.append(z)
                w[S] = ws
                Hs = [T] + [E.new() for j in range(p)]
                for u in range(N):
                    if not (S >> u) & 1: continue
                    G = unary(E, [AND(ws[x], both[x][u]) for x in range(N)], min(N, p + 1))
                    Gs = [T] + G
                    Hp = H[S ^ (1 << u)]
                    for j in range(p + 1):
                        for i in range(len(Gs)):
                            if i == 0: continue
                            if j + i >= p + 1: E.add([-Hp[j], -Gs[i]])
                            else: E.add([-Hp[j], -Gs[i], Hs[j + i]])
                    for j in range(1, p + 1): E.add([-Hp[j], Hs[j]])   # i = 0 transition
                H[S] = Hs
    PEEL = int(os.environ.get('PEEL', '0'))
    if PEEL > 0:
        # [P] iterated peeling: dim Fix(Theta) <= |core|, core = cells surviving repeated removal of cells with no in-edge or
        # no out-edge (support level).  alive_t literals imply true survival; require |alive_T| >= N^2 - p.
        n2 = N * N
        edge = {}
        for a in range(N):
            for b in range(N):
                for x in range(N):
                    for y in range(N):
                        lits = [AND(e[l][x][b], c[l][y][a]) for l in range(q)]
                        z = E.new(); E.add([-z] + lits)
                        for t in lits: E.add([z, -t])
                        edge[(a * N + b, x * N + y)] = z      # cell (a,b) -> cell (x,y)
        T_ = E.new(); E.add([T_])
        alive = [T_] * n2
        for t in range(PEEL):
            nxt = []
            for cc in range(n2):
                z = E.new(); E.add([-z, alive[cc]])
                ins = []; outs_ = []
                for d in range(n2):
                    if alive[d] is None: continue
                    w1 = E.new(); E.add([-w1, edge[(d, cc)]]); E.add([-w1, alive[d]]); ins.append(w1)
                    w2 = E.new(); E.add([-w2, edge[(cc, d)]]); E.add([-w2, alive[d]]); outs_.append(w2)
                E.add([-z] + ins); E.add([-z] + outs_)
                nxt.append(z)
            alive = nxt
        need = N * N - p
        cc_ = CardEnc.atleast(lits=alive, bound=need, top_id=E.nv, encoding=EncType.seqcounter)
        E.nv = max(E.nv, cc_.nv)
        for x in cc_.clauses: E.add(x)
    if extra is not None and r is not None and os.environ.get('TOP1') == '1':
        # term 1 = the term with |E_l| = r maximising |C_l| (valid: terms with equal |E| are interchangeable)
        s1 = sum(map(sum, extra))
        if s1 < N:
            for l in range(1, q):
                oE = unary(E, [e[l][i][j] for i in range(N) for j in range(N)], N)
                oC = unary(E, [c[l][i][j] for i in range(N) for j in range(N)], N)
                E.add([-oE[r - 1], -oC[s1]])
    if os.environ.get('NONZERO') == '1':   # every term nonzero (a zero term reduces to q-1, decided separately)
        for l in range(q):
            E.add([e[l][i][j] for i in range(N) for j in range(N)]); E.add([c[l][i][j] for i in range(N) for j in range(N)])
    if os.environ.get('TRANS') == '1':   # ansatz C_l = E_l^T (supports), fixed set of flip o transpose
        for l in range(q):
            for y in range(N):
                for v in range(N): E.add([-c[l][y][v], e[l][v][y]]); E.add([c[l][y][v], -e[l][v][y]])
    AUT = os.environ.get('AUT')
    if AUT:   # ansatz: (A,B,sigma) support automorphism: E_sig(l)[A x][B u] = E_l[x][u], C_sig(l)[B y][A v] = C_l[y][v]
        Ap, Bp, sg = [[int(t) for t in part.split(',')] for part in AUT.split(';')]
        def eq(a_, b_): E.add([-a_, b_]); E.add([a_, -b_])
        FL = os.environ.get('AUTF') == '1'   # flip type: E_sig(l)[x][u] = C_l[A x][B u], C_sig(l)[y][v] = E_l[B y][A v]
        for l in range(q):
            for x in range(N):
                for u in range(N):
                    if FL: eq(e[sg[l]][x][u], c[l][Ap[x]][Bp[u]]); eq(c[sg[l]][x][u], e[l][Bp[x]][Ap[u]])
                    else: eq(e[sg[l]][Ap[x]][Bp[u]], e[l][x][u]); eq(c[sg[l]][Bp[x]][Ap[u]], c[l][x][u])
    if extra is not None:   # fix C_1 = extra (0/1 matrix)
        for y in range(N):
            for v in range(N): E.add([c[0][y][v]] if extra[y][v] else [-c[0][y][v]])
    if sym:
        z, PV = rref(E, V, p, M)
        if LINK:   # U is determined by M at the pivot columns: U[a][k] = M[a][piv_k]; zero V rows -> U col 0
            for k in range(p):
                for a in range(M):
                    E.add([-z[k][M], -U[a][k]])
                    for j in range(M):
                        E.add([-PV[k][j], -U[a][k], Mv[(a, j)]]); E.add([-PV[k][j], U[a][k], -Mv[(a, j)]])
    return E, e, c, U, V

def xor_cnf(E):
    cl = list(E.cl); nv = E.nv
    for terms, rhs in E.xors:
        lits = terms[:]
        while len(lits) > 4:
            a, b, cc = lits[0], lits[1], lits[2]
            nv += 1; t = nv
            for bits in itertools.product([0, 1], repeat=3):
                par = sum(bits) % 2
                cl.append([(-l if bt else l) for l, bt in zip((a, b, cc), bits)] + ([t] if par else [-t]))
            lits = [t] + lits[3:]
        for bits in itertools.product([0, 1], repeat=len(lits)):
            if sum(bits) % 2 != rhs:
                cl.append([(-l if bt else l) for l, bt in zip(lits, bits)])
    return nv, cl

# ---------- exact Q rank ----------
def Kmat(N, Es, Cs):
    K = np.zeros((N * N, N * N), dtype=np.int64)
    for Em, Cm in zip(Es, Cs): K += np.kron(np.array(Em, dtype=np.int64), np.array(Cm, dtype=np.int64))
    return K
def swap(N):
    S = np.zeros((N * N, N * N), dtype=np.int64)
    for x in range(N):
        for y in range(N): S[x * N + y, y * N + x] = 1
    return S
def rankQ(N, Es, Cs):
    import flint
    Mx = swap(N) - Kmat(N, Es, Cs)
    return flint.fmpz_mat(Mx.tolist()).rank()
def rank_mod(Mx, P):
    A = (Mx % P).astype(np.int64); r = 0; rows, cols = A.shape
    for col in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i, col]: piv = i; break
        if piv is None: continue
        A[[r, piv]] = A[[piv, r]]; inv = pow(int(A[r, col]), P - 2, P); A[r] = (A[r] * inv) % P
        for i in range(rows):
            if i != r and A[i, col]: A[i] = (A[i] - A[i, col] * A[r]) % P
        r += 1
    return r

# ---------- stage 2: signs over GF(3) for a fixed support ----------
def sign_search(N, p, Es, Cs, maxiter=200, tl=None):
    """Es, Cs: 0/1 support matrices. Returns (signed Es, Cs, rankQ) with rankQ <= p, or None if no GF(3) signing."""
    from pysat.solvers import Solver
    q = len(Es); M = N * N
    # stage A: local search on signs (float SVD objective), exact Q check
    rng = np.random.default_rng(0); Sf = swap(N).astype(float)
    pos = [(0, l, x, u) for l in range(q) for x in range(N) for u in range(N) if Es[l][x][u]] + \
          [(1, l, x, u) for l in range(q) for x in range(N) for u in range(N) if Cs[l][x][u]]
    def bld(sg):
        Ef = np.array(Es, dtype=float); Cf = np.array(Cs, dtype=float)
        for s_, (side, l, x, u) in zip(sg, pos): (Ef if side == 0 else Cf)[l, x, u] = s_
        return Ef, Cf
    def fobj(sg):
        Ef, Cf = bld(sg); Kf = sum(np.kron(Ef[l], Cf[l]) for l in range(q))
        return np.linalg.svd(Sf - Kf, compute_uv=False)[p:].sum()
    tA = time.time()
    while time.time() - tA < float(os.environ.get('LSTIME', '15')) and pos:
        sg = rng.choice([-1., 1.], len(pos)); cur = fobj(sg)
        for it in range(2000):
            i = rng.integers(len(pos)); sg[i] *= -1; v = fobj(sg)
            if v <= cur + 1e-12: cur = v
            else: sg[i] *= -1
            if cur < 1e-9: break
        if cur < 1e-9:
            Ef, Cf = bld(sg); E2 = Ef.astype(int).tolist(); C2 = Cf.astype(int).tolist(); rq = rankQ(N, E2, C2)
            if rq <= p: return E2, C2, rq
    pool = [0]
    def new(): pool[0] += 1; return pool[0]
    cl = []
    sE = [[[new() if Es[l][x][u] else None for u in range(N)] for x in range(N)] for l in range(q)]
    sC = [[[new() if Cs[l][y][v] else None for v in range(N)] for y in range(N)] for l in range(q)]
    # gauge fixing: (A,B) diagonal signs + term signs: fix sign of a spanning forest... keep it simple: none (SAT side)
    def tri():
        a = [new(), new(), new()]; cl.append(a[:])
        for i, j in itertools.combinations(range(3), 2): cl.append([-a[i], -a[j]])
        return a
    U = [[tri() for k in range(p)] for a in range(M)]
    V = [[tri() for b in range(M)] for k in range(p)]
    def mul(a, b):
        c = tri()
        for i in range(3):
            for j in range(3): cl.append([-a[i], -b[j], c[(i * j) % 3]])
        return c
    def add(a, b):
        c = tri()
        for i in range(3):
            for j in range(3): cl.append([-a[i], -b[j], c[(i + j) % 3]])
        return c
    def pm(s):  # GF3 element from sign var: value 1 if not s, 2 if s
        c = tri(); cl.append([-c[0]]); cl.append([s, c[1]]); cl.append([-s, c[2]]); return c
    def xor(a, b):
        z = new(); cl.extend([[-z, a, b], [-z, -a, -b], [z, -a, b], [z, a, -b]]); return z
    Sw = swap(N)
    for x, y, u, v in itertools.product(range(N), repeat=4):
        terms = [mul(U[x * N + y][k], V[k][u * N + v]) for k in range(p)]
        for l in range(q):
            if sE[l][x][u] and sC[l][y][v]: terms.append(pm(xor(sE[l][x][u], sC[l][y][v])))
        acc = terms[0]
        for t in terms[1:]: acc = add(acc, t)
        cl.append([acc[int(Sw[x * N + y, u * N + v])]])
    S = Solver(name='cadical153', bootstrap_with=cl)
    signvars = [v for l in range(q) for X in (sE[l], sC[l]) for row in X for v in row if v]
    for it in range(maxiter):
        if not S.solve(): return None
        ms = set(l for l in S.get_model() if l > 0)
        sg = lambda v: -1 if v in ms else 1
        E2 = [[[sg(sE[l][x][u]) if sE[l][x][u] else 0 for u in range(N)] for x in range(N)] for l in range(q)]
        C2 = [[[sg(sC[l][y][v]) if sC[l][y][v] else 0 for v in range(N)] for y in range(N)] for l in range(q)]
        rq = rankQ(N, E2, C2)
        if rq <= p: return E2, C2, rq
        S.add_clause([(-v if v in ms else v) for v in signvars])
    return 'MAXITER'

UNSOUND = [0]
def run(N, p, q, r, tl=None, order=False, maxcegar=100000, out=None):
    import pycryptosat
    E, e, c, U, V = build(N, p, q, r, order=order)
    S = pycryptosat.Solver(threads=1, time_limit=tl) if tl else pycryptosat.Solver(threads=1)
    for cl in E.cl: S.add_clause(cl)
    for t, rr in E.xors: S.add_xor_clause(t, bool(rr))
    ev = [v for l in range(q) for X in (e[l], c[l]) for row in X for v in row]
    t0 = time.time(); nblock = 0
    while True:
        if tl:
            rem = tl - (time.time() - t0)
            if rem <= 0: return 'TIMEOUT', nblock, None
        sat, sol = S.solve()
        if sat is None: return 'TIMEOUT', nblock, None
        if not sat: return 'UNSAT', nblock, None
        val = lambda i: 1 if sol[i] else 0
        Es = [[[val(e[l][x][u]) for u in range(N)] for x in range(N)] for l in range(q)]
        Cs = [[[val(c[l][y][v]) for v in range(N)] for y in range(N)] for l in range(q)]
        res = sign_search(N, p, Es, Cs)
        if res is not None and res != 'MAXITER':
            E2, C2, rq = res
            return 'SAT', nblock, {"N": N, "p": p, "q": q, "rankQ": rq, "E": E2, "C": C2}
        nblock += 1
        if res == 'MAXITER': UNSOUND[0] += 1
        print(f"   support {nblock} unsignable ({res}) [{time.time()-t0:.1f}s]", flush=True)
        S.add_clause([(-v if sol[v] else v) for v in ev])
        if nblock >= maxcegar: return 'CEGAR-LIMIT', nblock, None

if __name__ == "__main__":
    N, p, q = map(int, sys.argv[1:4])
    rs = sys.argv[4] if len(sys.argv) > 4 else 'all'
    tl = float(os.environ.get('TL', '0')) or None
    order = bool(int(os.environ.get('ORD', '0')))
    rlist = range(N, 0, -1) if rs == 'all' else ([None] if rs == 'none' else [int(t) for t in rs.split(',')])
    T0 = time.time(); allunsat = True
    for r in rlist:
        t0 = time.time()
        st, nb, dec = run(N, p, q, r, tl=tl, order=order)
        print(f"MC N={N} (p,q)=({p},{q}) r={r}: {st} (blocked supports {nb}) [{time.time()-t0:.1f}s]", flush=True)
        if st == 'SAT':
            fn = f"mc_N{N}_{p}_{q}.json"; json.dump(dec, open(fn, 'w'))
            print("  rankQ =", dec['rankQ'], "saved", fn, flush=True); allunsat = False; break
        if st != 'UNSAT': allunsat = False
    if UNSOUND[0]: print("WARNING: %d supports blocked after MAXITER sign attempts -> UNSAT NOT a proof" % UNSOUND[0]); allunsat = False
    if allunsat: print(f"MC N={N} (p,q)=({p},{q}): UNSAT over GF(2) for all r => infeasible in MC [{time.time()-T0:.1f}s]")
