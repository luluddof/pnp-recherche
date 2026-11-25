# Cyclic ansatz (SAT direction only; negatives mean nothing): E_l = A^l E_0 B^-l, C_l = B^l C_0 A^-l (supports),
# A, B permutations with A^q = B^q = 1 (term l -> l+1 is an automorphism). GF(2) UV-SAT on (E_0, C_0) supports,
# then free GF(3) sign search on the resulting support and exact Q rank.
# usage: python3 ans.py N p q "cyclesA" "cyclesB"   e.g. "0,1,2;3,4,5"  ("" = identity)
import sys, os, itertools, time, json
import mc

def perm_from(cyc, N):
    P = list(range(N))
    if cyc:
        for cy in cyc.split(';'):
            c = [int(t) for t in cy.split(',')]
            for i in range(len(c)): P[c[i]] = c[(i + 1) % len(c)]
    return P

def powp(P, k):
    N = len(P); Q = list(range(N))
    for _ in range(k): Q = [P[Q[i]] for i in range(N)]
    return Q

def main():
    N, p, q = map(int, sys.argv[1:4]); A = perm_from(sys.argv[4], N); B = perm_from(sys.argv[5], N)
    assert powp(A, q) == list(range(N)) and powp(B, q) == list(range(N))
    TL = float(os.environ.get('TL', '540'))
    E = mc.Enc(); M = N * N
    e0 = [[E.new() for u in range(N)] for x in range(N)]
    c0 = [[E.new() for v in range(N)] for y in range(N)]
    for X in (e0, c0):
        for i in range(N):
            mc.amo(E, [X[i][j] for j in range(N)]); mc.amo(E, [X[j][i] for j in range(N)])
    E.add([e0[i][j] for i in range(N) for j in range(N)]); E.add([c0[i][j] for i in range(N) for j in range(N)])
    # E_l[A^l x][B^l u] = E_0[x][u] ; C_l[B^l y][A^l v] = C_0[y][v]
    e = []; c = []
    for l in range(q):
        Al, Bl = powp(A, l), powp(B, l)
        el = [[None] * N for _ in range(N)]; cl = [[None] * N for _ in range(N)]
        for x in range(N):
            for u in range(N):
                el[Al[x]][Bl[u]] = e0[x][u]; cl[Bl[x]][Al[u]] = c0[x][u]
        e.append(el); c.append(cl)
    U = [[E.new() for k in range(p)] for a in range(M)]
    V = [[E.new() for b in range(M)] for k in range(p)]
    def AND(a, b):
        z = E.new(); E.add([-z, a]); E.add([-z, b]); E.add([z, -a, -b]); return z
    for x, y, u, v in itertools.product(range(N), repeat=4):
        terms = [AND(U[x * N + y][k], V[k][u * N + v]) for k in range(p)]
        terms += [AND(e[l][x][u], c[l][y][v]) for l in range(q)]
        E.xors.append((terms, int(x == v and y == u)))
    mc.rref(E, V, p, M)
    import pycryptosat
    S = pycryptosat.Solver(threads=1, time_limit=TL)
    for x in E.cl: S.add_clause(x)
    for t, rr in E.xors: S.add_xor_clause(t, bool(rr))
    ev = [v for X in (e0, c0) for row in X for v in row]
    t0 = time.time(); nb = 0
    while True:
        sat, sol = S.solve()
        if sat is None: print(f"ANSATZ N={N} ({p},{q}) A={sys.argv[4]} B={sys.argv[5]}: TIMEOUT [{time.time()-t0:.1f}s]"); return
        if not sat: print(f"ANSATZ N={N} ({p},{q}) A={sys.argv[4]} B={sys.argv[5]}: UNSAT (ansatz only) blocked {nb} [{time.time()-t0:.1f}s]"); return
        Es = [[[int(sol[e[l][x][u]]) for u in range(N)] for x in range(N)] for l in range(q)]
        Cs = [[[int(sol[c[l][y][v]]) for v in range(N)] for y in range(N)] for l in range(q)]
        res = mc.sign_search(N, p, Es, Cs)
        if res is not None and res != 'MAXITER':
            E2, C2, rq = res
            fn = f"ans_N{N}_{p}_{q}.json"; json.dump({"N": N, "p": p, "q": q, "rankQ": rq, "E": E2, "C": C2}, open(fn, 'w'))
            print(f"ANSATZ N={N} ({p},{q}) A={sys.argv[4]} B={sys.argv[5]}: SAT rankQ {rq} saved {fn} [{time.time()-t0:.1f}s]"); return
        nb += 1; print(f"  support {nb} unsignable ({res})", flush=True)
        S.add_clause([(-v if sol[v] else v) for v in ev])

if __name__ == "__main__":
    main()
