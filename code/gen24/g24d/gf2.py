# G24d: exact feasibility of (p,q) over GF(2) with GL_p x GL_q symmetry breaking (RREF).
# R = U V (V: p x N^2 in reduced row echelon form, zero rows last) ; C = Gm Hm (Hm: q x N^2 in RREF, (y,v)-side).
# The RREF of a row space is unique, so this removes the GL_p(F2) x GL_q(F2) factor-basis symmetry exactly.
# Optional: XOR mode (CryptoMiniSat native xor) or pure CNF (Tseitin xor chains) for CaDiCaL/kissat + DRAT.
import sys, itertools, time, json, os

class Enc:
    def __init__(self): self.nv = 0; self.cl = []; self.xors = []
    def new(self): self.nv += 1; return self.nv
    def add(self, c): self.cl.append(c)

def rref(E, A, k, m):
    """A[i][j] vars, k rows, m cols: constrain to RREF with zero rows last."""
    z = [[None] * (m + 1) for i in range(k)]
    T = E.new(); E.add([T])
    for i in range(k):
        z[i][0] = T
        for j in range(m):
            v = E.new(); z[i][j + 1] = v       # v <-> z[i][j] & -A[i][j]
            E.add([-v, z[i][j]]); E.add([-v, -A[i][j]]); E.add([v, -z[i][j], A[i][j]])
    for i in range(k):
        for j in range(m):
            pv = E.new()                       # pv <-> z[i][j] & A[i][j]
            E.add([-pv, z[i][j]]); E.add([-pv, A[i][j]]); E.add([pv, -z[i][j], -A[i][j]])
            if i + 1 < k: E.add([-pv, z[i + 1][j + 1]])
            for i2 in range(k):
                if i2 != i: E.add([-pv, -A[i2][j]])
        if i + 1 < k: E.add([-z[i][m], z[i + 1][m]])

def build(N, p, q, sym=True):
    E = Enc(); M = N * N
    U = [[E.new() for k in range(p)] for a in range(M)]
    V = [[E.new() for b in range(M)] for k in range(p)]
    G = [[E.new() for l in range(q)] for a in range(M)]   # G[(x,u)][l]
    H = [[E.new() for b in range(M)] for l in range(q)]   # H[l][(y,v)]
    prod = {}
    def AND(a, b):
        key = (min(a, b), max(a, b))
        if key not in prod:
            z = E.new(); prod[key] = z
            E.add([-z, a]); E.add([-z, b]); E.add([z, -a, -b])
        return prod[key]
    for x, y, u, v in itertools.product(range(N), repeat=4):
        terms = [AND(U[x * N + y][k], V[k][u * N + v]) for k in range(p)]
        terms += [AND(G[x * N + u][l], H[l][y * N + v]) for l in range(q)]
        E.xors.append((terms, int(x == v and y == u)))
    r1 = os.environ.get('GNF')
    if r1 is not None:   # G_1 = diag(I_r1, 0) via GL_N x GL_N (x,v by A,A^-T; y,u by B,B^-T); RREF on V only
        r1 = int(r1)
        for x in range(N):
            for u in range(N): E.add([G[x * N + u][0]] if (x == u and x < r1) else [-G[x * N + u][0]])
        if sym: rref(E, V, p, M)
    elif sym:
        rref(E, V, p, M); rref(E, H, q, M)
    return E, U, V, G, H

def xor_cnf(E):
    """Tseitin-chain the xors into CNF (cut length 4)."""
    cl = list(E.cl)
    nv = E.nv
    for terms, rhs in E.xors:
        lits = terms[:]
        while len(lits) > 4:
            a, b, c = lits[0], lits[1], lits[2]
            nv += 1; t = nv     # t = a ^ b ^ c
            for bits in itertools.product([0, 1], repeat=3):
                par = sum(bits) % 2
                # forbid assignments where t != parity
                cl.append([(-l if bt else l) for l, bt in zip((a, b, c), bits)] + ([t] if par else [-t]))
            lits = [t] + lits[3:]
        k = len(lits)
        for bits in itertools.product([0, 1], repeat=k):
            if sum(bits) % 2 != rhs:
                cl.append([(-l if bt else l) for l, bt in zip(lits, bits)])
    return nv, cl

def solve(N, p, q, sym=True, mode='cms', dimacs=None, tl=None):
    E, U, V, G, H = build(N, p, q, sym)
    if mode == 'cms':
        import pycryptosat
        S = pycryptosat.Solver(threads=1, time_limit=tl) if tl else pycryptosat.Solver(threads=1)
        for c in E.cl: S.add_clause(c)
        for t, r in E.xors: S.add_xor_clause(t, bool(r))
        sat, sol = S.solve()
        if sat is None: return 'UNKNOWN', None
        if not sat: return 'UNSAT', None
        val = lambda i: 1 if sol[i] else 0
    else:
        nv, cl = xor_cnf(E)
        if dimacs:
            with open(dimacs, 'w') as f:
                f.write(f"p cnf {nv} {len(cl)}\n")
                for c in cl: f.write(' '.join(map(str, c)) + ' 0\n')
            if mode == 'dimacs': return 'WRITTEN', None
        from pysat.solvers import Solver
        S = Solver(name=mode, bootstrap_with=cl)
        if not S.solve(): return 'UNSAT', None
        ms = set(l for l in S.get_model() if l > 0)
        val = lambda i: 1 if i in ms else 0
    dec = {"N": N, "p": p, "q": q, "field": 2,
           "U": [[val(U[a][k]) for k in range(p)] for a in range(N * N)],
           "V": [[val(V[k][b]) for b in range(N * N)] for k in range(p)],
           "G": [[[val(G[x * N + u][l]) for u in range(N)] for x in range(N)] for l in range(q)],
           "H": [[[val(H[l][y * N + v]) for v in range(N)] for y in range(N)] for l in range(q)]}
    return 'SAT', dec

if __name__ == "__main__":
    N, p, q = map(int, sys.argv[1:4]); mode = sys.argv[4] if len(sys.argv) > 4 else 'cms'
    sym = not os.environ.get('NOSYM')
    t0 = time.time(); st, dec = solve(N, p, q, sym, mode, dimacs=os.environ.get('DIMACS'))
    print(f"GF(2) N={N} (p,q)=({p},{q}) mode={mode} sym={sym}: {st} [{time.time()-t0:.1f}s]", flush=True)
    if dec:
        fn = f"dec_gf2_N{N}_{p}_{q}.json"; json.dump(dec, open(fn, 'w'))
        import verify_gf
        print("  verify_gf: rr(R)=%d rc(C)=%d (re-expansion OK)" % verify_gf.verify(dec), flush=True)
