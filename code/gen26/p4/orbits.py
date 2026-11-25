# orbit representatives of C_1 (0/1 partial perm, rank <= r) under the stabiliser of E_1 = diag(I_r,0):
# (y,v) -> (B y, A v), A = pi+alpha, B = pi+beta (pi on [r], alpha, beta on complement), plus transpose (y,v)->(v,y).
import itertools
def reps(N, r):
    co = list(range(r, N))
    G = []
    for pi in itertools.permutations(range(r)):
        for al in itertools.permutations(co):
            for be in itertools.permutations(co):
                A = list(pi) + list(al); B = list(pi) + list(be); G.append((A, B))
    seen = set(); out = []
    for k in range(0, r + 1):
        for rows in itertools.combinations(range(N), k):
            for cols in itertools.permutations(range(N), k):
                C = frozenset(zip(rows, cols))
                if C in seen: continue
                orb = set()
                for A, B in G:
                    D = frozenset((B[y], A[v]) for y, v in C); orb.add(D); orb.add(frozenset((v, y) for y, v in D))
                seen |= orb
                out.append((len(orb), sorted(min(orb, key=lambda s: sorted(s)))))
    return out
if __name__ == "__main__":
    import sys
    N, r = map(int, sys.argv[1:3]); R = reps(N, r); print(N, r, len(R), "orbits; total", sum(o for o, _ in R))
