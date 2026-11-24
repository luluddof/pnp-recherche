# Exact feasibility of (p,q) for the crossing over GF(P), P an odd prime (or 2), by CNF (one-hot field elements, mod-P
# running sums) + CaDiCaL.  Same unknowns as sat2.py.  SAT answers are re-verified by verify_gf.py.
import sys, itertools, time, json
from pysat.solvers import Cadical153
from pysat.card import CardEnc
def solve(N, p, q, P, budget=None):
    cl = []; nv = [0]
    def new():
        nv[0] += 1; return nv[0]
    def elem():
        e = [new() for _ in range(P)]
        cl.extend(CardEnc.equals(lits=e, bound=1, top_id=10**7).clauses if False else exactly_one(e))
        return e
    def exactly_one(e):
        c = [list(e)]
        for i in range(len(e)):
            for j in range(i+1, len(e)): c.append([-e[i], -e[j]])
        return c
    U = {(a, k): elem() for a in range(N*N) for k in range(p)}
    V = {(k, b): elem() for k in range(p) for b in range(N*N)}
    G = {(l, x, u): elem() for l in range(q) for x in range(N) for u in range(N)}
    H = {(l, y, v): elem() for l in range(q) for y in range(N) for v in range(N)}
    pc = {}
    def MUL(a, b):  # returns one-hot list of P literals for a*b
        key = (id(a), id(b))
        if key in pc: return pc[key]
        c = [new() for _ in range(P)]; cl.extend(exactly_one(c))
        for i in range(P):
            for j in range(P):
                cl.append([-a[i], -b[j], c[(i*j) % P]])
        pc[key] = c; return c
    for x, y, u, v in itertools.product(range(N), repeat=4):
        terms = [MUL(U[(x*N+y, k)], V[(k, u*N+v)]) for k in range(p)] + [MUL(G[(l, x, u)], H[(l, y, v)]) for l in range(q)]
        rhs = 1 if (x == v and y == u) else 0
        s = terms[0]
        for t in terms[1:]:
            s2 = [new() for _ in range(P)]; cl.extend(exactly_one(s2))
            for i in range(P):
                for j in range(P): cl.append([-s[i], -t[j], s2[(i+j) % P]])
            s = s2
        cl.append([s[rhs]])
    # nonzero terms (monotone in p, q)
    for l in range(q):
        cl.append([-G[(l, x, u)][0] for x in range(N) for u in range(N)])
        cl.append([-H[(l, y, v)][0] for y in range(N) for v in range(N)])
    for k in range(p):
        cl.append([-U[(a, k)][0] for a in range(N*N)]); cl.append([-V[(k, b)][0] for b in range(N*N)])
    S = Cadical153(bootstrap_with=cl)
    ok = S.solve()
    if not ok: return "UNSAT", None
    m = set(l for l in S.get_model() if l > 0)
    val = lambda e: next(i for i in range(P) if e[i] in m)
    dec = {"N": N, "p": p, "q": q, "field": P,
           "U": [[val(U[(a, k)]) for k in range(p)] for a in range(N*N)],
           "V": [[val(V[(k, b)]) for b in range(N*N)] for k in range(p)],
           "G": [[[val(G[(l, x, u)]) for u in range(N)] for x in range(N)] for l in range(q)],
           "H": [[[val(H[(l, y, v)]) for v in range(N)] for y in range(N)] for l in range(q)]}
    return "SAT", dec
if __name__ == "__main__":
    N, p, q, P = map(int, sys.argv[1:5])
    t0 = time.time(); st, dec = solve(N, p, q, P)
    print(f"GF({P}) N={N} (p,q)=({p},{q}): {st}  [{time.time()-t0:.1f}s]", flush=True)
    if dec: json.dump(dec, open(f"dec_gf{P}_N{N}_{p}_{q}.json", "w"))
