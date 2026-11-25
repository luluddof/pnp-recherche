# support-only feasibility (implied constraints only: partial perms, normal form, RSUM, CORE, STAIR) - no rank encoding
import sys, os, time
os.environ['NOUV'] = '1'
import mc, pycryptosat
N, p, q, r = map(int, sys.argv[1:5])
E, e, c, U, V = mc.build(N, p, q, r, order=True)
S = pycryptosat.Solver(threads=1, time_limit=200)
for x in E.cl: S.add_clause(x)
t0 = time.time(); sat, sol = S.solve()
print(f"support-only N={N} p={p} q={q} r={r}: {sat} [{time.time()-t0:.1f}s] vars {E.nv} clauses {len(E.cl)}")
