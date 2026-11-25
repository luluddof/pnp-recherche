import sys, os, time
os.environ['NOUV'] = '1'
import mc, pycryptosat
N, p, q, r = map(int, sys.argv[1:5]); TL = float(sys.argv[5])
E, e, c, U, V = mc.build(N, p, q, r, order=True)
S = pycryptosat.Solver(threads=1)
for x in E.cl: S.add_clause(x)
ev = [v for l in range(q) for X in (e[l], c[l]) for row in X for v in row]
t0 = time.time(); n = 0
while time.time() - t0 < TL:
    sat, sol = S.solve()
    if not sat: print("exhausted", n); break
    n += 1; S.add_clause([(-v if sol[v] else v) for v in ev])
print(f"N={N} p={p} q={q} r={r}: {n} support solutions in {time.time()-t0:.1f}s")
