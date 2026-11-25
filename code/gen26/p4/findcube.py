# which C_1-orbit cube contains a known solution (json with E, C lists)
import sys, json, numpy as np, orbits, itertools
fn, keyE, keyC = sys.argv[1], sys.argv[2], sys.argv[3]
d = json.load(open(fn)); Es = [np.abs(np.array(m, dtype=int)) for m in d[keyE]]; Cs = [np.abs(np.array(m, dtype=int)) for m in d[keyC]]
N = Es[0].shape[0]; r = max(int(m.sum()) for m in Es + Cs)
R = orbits.reps(N, r); idx = {}
co = list(range(r, N)); G = []
for pi in itertools.permutations(range(r)):
    for al in itertools.permutations(co):
        for be in itertools.permutations(co): G.append((list(pi) + list(al), list(pi) + list(be)))
for i, (o, rep) in enumerate(R):
    C = frozenset(map(tuple, rep))
    for A, B in G:
        D = frozenset((B[y], A[v]) for y, v in C); idx[D] = i; idx[frozenset((v, y) for y, v in D)] = i
for flip in (0, 1):
    for l in range(len(Es)):
        E1, C1 = (Es[l], Cs[l]) if not flip else (Cs[l], Es[l])
        if E1.sum() != r: continue
        pairs = [(x, u) for x in range(N) for u in range(N) if E1[x, u]]
        rows = [x for x, u in pairs] + [x for x in range(N) if x not in [a for a, b in pairs]]
        cols = [u for x, u in pairs] + [u for u in range(N) if u not in [b for a, b in pairs]]
        # A maps row rows[i] -> i ; B maps col cols[i] -> i ; E1 -> A E1 B^T = diag
        A = {rows[i]: i for i in range(N)}; B = {cols[i]: i for i in range(N)}
        # C1 -> B C1 A^T : entry (y,v) -> (B y, A v)
        D = frozenset((B[y], A[v]) for y in range(N) for v in range(N) if C1[y, v])
        print("flip", flip, "term", l, "r", r, "cube", idx[D], R[idx[D]][1])
