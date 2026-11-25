# sample support-only solutions (all lemmas) and report their GF(2) ranks
import sys, os, time, numpy as np
os.environ['NOUV'] = '1'
import mc, pycryptosat
N, p, q, r, K = map(int, sys.argv[1:6])
E, e, c, U, V = mc.build(N, p, q, r, order=True)
S = pycryptosat.Solver(threads=1)
for x in E.cl: S.add_clause(x)
ev = [v for l in range(q) for X in (e[l], c[l]) for row in X for v in row]
t0 = time.time(); ranks = []
for k in range(K):
    sat, sol = S.solve()
    if not sat: print("exhausted"); break
    Es = [[[int(sol[e[l][x][u]]) for u in range(N)] for x in range(N)] for l in range(q)]
    Cs = [[[int(sol[c[l][y][v]]) for v in range(N)] for y in range(N)] for l in range(q)]
    ranks.append(mc.rank_mod((mc.swap(N) + mc.Kmat(N, Es, Cs)) % 2, 2))
    S.add_clause([(-v if sol[v] else v) for v in ev])
print(f"N={N} p={p} q={q} r={r}: {len(ranks)} samples, GF(2) ranks min {min(ranks)} mean {np.mean(ranks):.1f} hist {np.bincount(ranks)[p:]} [{time.time()-t0:.0f}s]")
