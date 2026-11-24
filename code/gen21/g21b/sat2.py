# Exact feasibility of (p,q) for the crossing X=[x=v][y=u] over GF(2) via CryptoMiniSat (native XOR clauses).
# Unknowns: R = U V  (U: N^2 x p rows (x,y), V: p x N^2 cols (u,v)),  C = sum_l G_l(x,u) H_l(y,v).
# Constraint (all N^4 entries):  sum_k U[xy,k]V[k,uv] + sum_l G_l[x,u]H_l[y,v] = X  (mod 2).
# SAT answer is re-verified by independent exact rank computation over GF(2) (verify_gf.py).
import sys, itertools, time, json, pycryptosat
def solve(N, p, q, seed=0, timeout=None, symbreak=True, r1=None):
    S = pycryptosat.Solver(threads=1, time_limit=timeout) if timeout else pycryptosat.Solver(threads=1)
    nv = [0]
    def new():
        nv[0] += 1; return nv[0]
    U = {(a, k): new() for a in range(N*N) for k in range(p)}
    V = {(k, b): new() for k in range(p) for b in range(N*N)}
    G = {(l, x, u): new() for l in range(q) for x in range(N) for u in range(N)}
    H = {(l, y, v): new() for l in range(q) for y in range(N) for v in range(N)}
    prod = {}
    def AND(a, b):
        key = (min(a, b), max(a, b))
        if key in prod: return prod[key]
        z = new(); prod[key] = z
        S.add_clause([-z, a]); S.add_clause([-z, b]); S.add_clause([z, -a, -b])
        return z
    for x, y, u, v in itertools.product(range(N), repeat=4):
        terms = [AND(U[(x*N+y, k)], V[(k, u*N+v)]) for k in range(p)]
        terms += [AND(G[(l, x, u)], H[(l, y, v)]) for l in range(q)]
        rhs = (x == v and y == u)
        if terms: S.add_xor_clause(terms, rhs)
        elif rhs: return "UNSAT", None
    if symbreak:
        # each G_l nonzero and each H_l nonzero (else drop the term; feasibility is monotone in q), same for U cols / V rows
        for l in range(q):
            S.add_clause([G[(l, x, u)] for x in range(N) for u in range(N)])
            S.add_clause([H[(l, y, v)] for y in range(N) for v in range(N)])
        for k in range(p):
            S.add_clause([U[(a, k)] for a in range(N*N)])
            S.add_clause([V[(k, b)] for b in range(N*N)])
    if r1 is not None:  # G_1 in rank normal form diag(I_r1, 0) (GL_N x GL_N symmetry of X: x,v by A,A^-1; y,u by B,B^-1)
        for x in range(N):
            for u in range(N):
                lit = G[(0, x, u)]
                S.add_clause([lit] if (x == u and x < r1) else [-lit])
    sat, sol = S.solve()
    if sat is None: return "UNKNOWN", None
    if not sat: return "UNSAT", None
    val = lambda i: 1 if sol[i] else 0
    dec = {"N": N, "p": p, "q": q, "field": 2,
           "U": [[val(U[(a, k)]) for k in range(p)] for a in range(N*N)],
           "V": [[val(V[(k, b)]) for b in range(N*N)] for k in range(p)],
           "G": [[[val(G[(l, x, u)]) for u in range(N)] for x in range(N)] for l in range(q)],
           "H": [[[val(H[(l, y, v)]) for v in range(N)] for y in range(N)] for l in range(q)]}
    return "SAT", dec
if __name__ == "__main__":
    N, p, q = map(int, sys.argv[1:4]); to = float(sys.argv[4]) if len(sys.argv) > 4 else None
    r1 = int(sys.argv[5]) if len(sys.argv) > 5 else None
    t0 = time.time(); st, dec = solve(N, p, q, timeout=to, r1=r1)
    print(f"GF(2) N={N} (p,q)=({p},{q}) r1={r1}: {st}  [{time.time()-t0:.1f}s]", flush=True)
    if dec: json.dump(dec, open(f"dec_gf2_N{N}_{p}_{q}.json", "w"))
