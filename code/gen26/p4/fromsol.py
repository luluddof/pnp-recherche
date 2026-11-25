# kissat/cadical model of a dumped GF(2) MC instance -> supports -> GF(3) sign search -> exact Q rank -> json
import sys, os, json, mc
N, p, q, r = map(int, sys.argv[1:5]); solf = sys.argv[5]
E, e, c, U, V = mc.build(N, p, q, r, order=bool(int(os.environ.get('ORD', '0'))))
ms = set()
for line in open(solf):
    if line.startswith('v'): ms.update(int(t) for t in line.split()[1:] if int(t) > 0)
Es = [[[int(e[l][x][u] in ms) for u in range(N)] for x in range(N)] for l in range(q)]
Cs = [[[int(c[l][y][v] in ms) for v in range(N)] for y in range(N)] for l in range(q)]
import numpy as np
M2 = (mc.swap(N) + mc.Kmat(N, Es, Cs)) % 2
print("GF(2) rank of support matrix:", mc.rank_mod(M2, 2))
res = mc.sign_search(N, p, Es, Cs)
if res is None or res == 'MAXITER': print("unsignable:", res); sys.exit(1)
E2, C2, rq = res
print("signed: rank_Q(Swap-K) =", rq, " rank mod 3:", mc.rank_mod(mc.swap(N) - mc.Kmat(N, E2, C2), 3))
fn = f"mc_N{N}_{p}_{q}.json"; json.dump({"N": N, "p": p, "q": q, "rankQ": rq, "E": E2, "C": C2}, open(fn, 'w')); print("saved", fn)
