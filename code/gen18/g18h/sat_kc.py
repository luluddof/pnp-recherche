"""G18h (c): GF(2), N = 4 (lambda = 2).  Is the crossing [x=v][y=u] a sum of k pieces with prescribed
(row rank <= r_i, col rank <= c_i)?  A piece (r, c): r row terms f(x,y)g(u,v)  AND the same function as c col
terms f'(x,u)g'(y,v).  r=None or c=None: that direction unconstrained.  Floor profile: all r_i*c_i <= N^2 = 16.
Model is re-verified exactly (numpy, GF(2) ranks) when SAT."""
import itertools, sys, time, threading
import numpy as np
from pysat.solvers import Cadical153, Glucose4

def gf2rank(M):
    M = M.copy() % 2; r = 0
    for c in range(M.shape[1]):
        p = None
        for i in range(r, M.shape[0]):
            if M[i, c]: p = i; break
        if p is None: continue
        M[[r, p]] = M[[p, r]]
        for i in range(M.shape[0]):
            if i != r and M[i, c]: M[i] ^= M[r]
        r += 1
    return r

def run(N, profile, tl, tag=""):
    nv = [0]
    def new():
        nv[0] += 1; return nv[0]
    cl = []
    idx = list(itertools.product(range(N), repeat=2))
    def AND(a, b):
        p = new(); cl.extend([[-p, a], [-p, b], [p, -a, -b]]); return p
    def XOR(lits, rhs):   # sum lits = rhs mod 2, chained
        lits = list(lits)
        if not lits:
            if rhs: cl.append([])
            return
        acc = lits[0]
        for l in lits[1:]:
            t = new()
            cl.extend([[-t, acc, l], [-t, -acc, -l], [t, -acc, l], [t, acc, -l]])
            acc = t
        cl.append([acc] if rhs else [-acc])
    pieces = []
    for (r, c) in profile:
        rowt = [({p: new() for p in idx}, {p: new() for p in idx}) for _ in range(r or 0)]
        colt = [({p: new() for p in idx}, {p: new() for p in idx}) for _ in range(c or 0)]
        pieces.append((r, c, rowt, colt))
    quads = list(itertools.product(range(N), repeat=4))
    for x, y, u, v in quads:
        tot = []
        for r, c, rowt, colt in pieces:
            rl = [AND(f[(x, y)], g[(u, v)]) for f, g in rowt]
            cll = [AND(f[(x, u)], g[(y, v)]) for f, g in colt]
            if r is not None and c is not None:
                XOR(rl + cll, 0)       # both factorizations give the same function
                tot += rl
            elif r is not None: tot += rl
            else: tot += cll
        XOR(tot, int(x == v and y == u))
    s = Cadical153(bootstrap_with=cl)
    t0 = time.time()
    timer = threading.Timer(tl, lambda: s.interrupt()); timer.start()
    try:
        res = s.solve()
    except Exception as e:
        res = None
    timer.cancel()
    out = {True: "SAT", False: "UNSAT", None: "timeout"}[res]
    print(f"{tag} N={N} profile={profile}: {out} ({time.time()-t0:.0f}s, vars {nv[0]}, clauses {len(cl)})"); sys.stdout.flush()
    if res:
        m = set(l for l in s.get_model() if l > 0)
        tot = np.zeros((N,) * 4, dtype=np.int64)
        for r, c, rowt, colt in pieces:
            A = np.zeros((N,) * 4, dtype=np.int64); B = np.zeros((N,) * 4, dtype=np.int64)
            for x, y, u, v in quads:
                A[x, y, u, v] = sum((f[(x, y)] in m) and (g[(u, v)] in m) for f, g in rowt) % 2
                B[x, y, u, v] = sum((f[(x, u)] in m) and (g[(y, v)] in m) for f, g in colt) % 2
            if r is not None and c is not None: assert (A == B).all()
            P = A if r is not None else B
            rr = gf2rank(P.reshape(N * N, N * N)); rc = gf2rank(P.transpose(0, 2, 1, 3).reshape(N * N, N * N))
            print(f"    piece: exact GF(2) ranks row {rr} col {rc}  product {rr*rc}")
            tot = (tot + P) % 2
        X = np.zeros((N,) * 4, dtype=np.int64)
        for x, y in idx: X[x, y, y, x] = 1
        print("    sum of pieces == crossing:", bool((tot == X).all()))
    s.delete(); return res

if __name__ == "__main__":
    tl = int(sys.argv[1])
    run(2, [(None, 1), (1, None)], 20, "sanity")
    run(4, [(None, 1), (1, None)], 60, "k=2")
    run(4, [(None, 1), (4, 4), (4, 4), (1, None)], tl, "k=4 (tensor square of the 1-bit identity exists)")
    for prof in sys.argv[2:]:
        run(4, eval(prof), tl, "k=3")
