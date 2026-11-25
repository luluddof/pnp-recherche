# CEGAR for MC(N,p,q) branch r over GF(2) WITHOUT the U V rank encoding.
# Variables: supports e,c (partial perms, normal form as mc.build), kt = AND(e,c) (edge literals), mv = Swap + K mod 2.
# Loop: solve; M = current GF(2) matrix; if rank_GF2(M) <= p -> sign stage (GF(3) SAT + exact Q rank) ;
#   else add valid lemmas:  (a) acyclic-set lemma: S (|S| = p+1 cells) with no backward/self edge of the cell digraph
#   in some order => (I+Theta)[S,S] unitriangular => rank > p.  Clause: OR of edge literals of backward/self pairs in S.
#   (b) minor lemma: R x C (p+1)x(p+1) nonsingular block of M => clause: some mv in the block changes.
# Both are implied by rank <= p, so UNSAT of the final formula => GF(2)-UNSAT => MC-UNSAT (for that branch).
import sys, os, time, itertools, json, random
import numpy as np
import mc

def gf2_rank(A):
    A = A.copy() % 2; r = 0; rows, cols = A.shape
    for c in range(cols):
        piv = np.nonzero(A[r:, c])[0]
        if len(piv) == 0: continue
        i = r + piv[0]; A[[r, i]] = A[[i, r]]
        nz = np.nonzero(A[:, c])[0]; nz = nz[nz != r]
        A[nz] ^= A[r]; r += 1
        if r == rows: break
    return r

def nonsing_minor(A, k, rng):
    """random greedy: rows/cols of a k x k nonsingular block of A over GF(2) (A has rank >= k)."""
    n, m = A.shape
    rows = list(rng.permutation(n)); cols = list(rng.permutation(m))
    B = A[rows][:, cols].copy() % 2
    R, C = [], []; r = 0
    # elimination with row pivoting over columns in random order; record original row/col of pivots
    rowid = list(rows)
    for jc in range(m):
        piv = np.nonzero(B[r:, jc])[0]
        if len(piv) == 0: continue
        i = r + piv[0]; B[[r, i]] = B[[i, r]]; rowid[r], rowid[i] = rowid[i], rowid[r]
        nz = np.nonzero(B[:, jc])[0]; nz = nz[nz > r]
        B[nz] ^= B[r]
        R.append(rowid[r]); C.append(cols[jc]); r += 1
        if r == k: break
    return R, C

def main():
    N, p, q, r = map(int, sys.argv[1:5]); TL = float(sys.argv[5]) if len(sys.argv) > 5 else 550
    os.environ['LINK'] = '1'
    E, e, c, U, V = mc.build(N, p, q, r, sym=False, order=bool(int(os.environ.get('ORD', '1'))), c22=True)
    # mc.build with LINK=1 created mv vars and xors (kt+[mv], swap) and (UV terms + [mv], 0); drop the UV xors
    import pycryptosat
    S = pycryptosat.Solver(threads=1)
    for cl in E.cl: S.add_clause(cl)
    M2 = N * N
    # recover kt and mv ids from the xor list: even entries = (kt + [mv], swap), odd = (UV + [mv], 0)
    kt = {}; mv = {}
    idx = 0
    for x, y, u, v in itertools.product(range(N), repeat=4):
        terms, rhs = E.xors[2 * idx]; idx += 1
        S.add_xor_clause(terms, bool(rhs))
        a, b = x * N + y, u * N + v
        mv[(a, b)] = terms[-1]; kt[(a, b)] = terms[:-1]
    # edge literal: Theta[(x,y),(a,b)] = K[(x,y),(b,a)] ; cell (a,b) -> (x,y) uses kt[(x*N+y, b*N+a)]
    def edge_lits(cfrom, cto):
        a, b = divmod(cfrom, N); x, y = divmod(cto, N)
        return kt[(x * N + y, b * N + a)]
    ev = [v for l in range(q) for X in (e[l], c[l]) for row in X for v in row]
    rng = np.random.default_rng(1); t0 = time.time(); it = 0; nlem = [0, 0]
    Sw = mc.swap(N)
    while True:
        if time.time() - t0 > TL: print(f"TIMEOUT after {it} iterations, lemmas acyc {nlem[0]} minor {nlem[1]}"); return 'TIMEOUT'
        sat, sol = S.solve()
        if not sat:
            print(f"UNSAT after {it} iterations, lemmas acyc {nlem[0]} minor {nlem[1]} [{time.time()-t0:.1f}s]"); return 'UNSAT'
        it += 1
        Es = [[[int(sol[e[l][x][u]]) for u in range(N)] for x in range(N)] for l in range(q)]
        Cs = [[[int(sol[c[l][y][v]]) for v in range(N)] for y in range(N)] for l in range(q)]
        M = (Sw + mc.Kmat(N, Es, Cs)) % 2
        rk = gf2_rank(M)
        if rk <= p:
            res = mc.sign_search(N, p, Es, Cs)
            if res is not None and res != 'MAXITER':
                E2, C2, rq = res
                fn = f"mc_N{N}_{p}_{q}.json"; json.dump({"N": N, "p": p, "q": q, "rankQ": rq, "E": E2, "C": C2}, open(fn, 'w'))
                print(f"SAT rankQ {rq} saved {fn} after {it} iterations [{time.time()-t0:.1f}s]"); return 'SAT'
            if res == 'MAXITER': print("SIGN MAXITER -> stop (not a proof)"); return 'ERR'
            S.add_clause([(-v if sol[v] else v) for v in ev]); continue   # provably unsignable support
        # Theta (GF(2)) adjacency of the current support: edge from cell cf to ct if any kt literal true
        Th = M.T @ np.zeros(1) if False else None
        out = [[] for _ in range(M2)]
        loop = [False] * M2
        for cf in range(M2):
            for ct in range(M2):
                if any(sol[l] for l in edge_lits(cf, ct)):
                    if cf == ct: loop[cf] = True
                    else: out[cf].append(ct)
        added = 0
        for trial in range(30):
            order = rng.permutation(M2); pos = np.empty(M2, int); pos[order] = np.arange(M2)
            # greedy acyclic set: scan cells in order, keep c if no edge between c and kept cells going backward
            Sset = []
            for cc in order:
                if loop[cc]: continue
                # edges c -> d with d earlier (pos[d] < pos[c]) inside S are backward; forbid them, and edges d -> c with
                # d in S are forward (allowed). since we scan in increasing pos, new c is last: need no edge c -> d in S
                if any(d in Sset for d in out[cc]): continue
                Sset.append(cc)
            if len(Sset) >= p + 1:
                Sset = Sset[:p + 1]
                lits = set()
                for i, ci in enumerate(Sset):
                    for j in range(i + 1):
                        lits.update(edge_lits(ci, Sset[j]))   # edges from later/equal cell to earlier/equal cell
                S.add_clause(sorted(lits)); added += 1; nlem[0] += 1
        if added == 0:
            for trial in range(10):
                R, C = nonsing_minor(M, p + 1, rng)
                S.add_clause([(-mv[(a, b)] if M[a, b] else mv[(a, b)]) for a in R for b in C]); nlem[1] += 1
        if it % 200 == 0: print(f"  it {it} rank {rk} lemmas acyc {nlem[0]} minor {nlem[1]} [{time.time()-t0:.1f}s]", flush=True)

if __name__ == "__main__":
    main()
