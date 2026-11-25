# extract supports from a model (e vars = first q*N*N vars, c vars next, as in every mc.build version), sign search, Q rank
import sys, json, mc
N, p, q = map(int, sys.argv[1:4]); solf = sys.argv[4]
ms = set()
for line in open(solf):
    if line.startswith('v'): ms.update(int(t) for t in line.split()[1:] if int(t) > 0)
idx = lambda l, x, u, side: 1 + side * q * N * N + l * N * N + x * N + u
Es = [[[int(idx(l, x, u, 0) in ms) for u in range(N)] for x in range(N)] for l in range(q)]
Cs = [[[int(idx(l, y, v, 1) in ms) for v in range(N)] for y in range(N)] for l in range(q)]
print("GF(2) rank:", mc.rank_mod((mc.swap(N) + mc.Kmat(N, Es, Cs)) % 2, 2))
res = mc.sign_search(N, p, Es, Cs)
if res is None or res == 'MAXITER': print("unsignable", res); sys.exit()
E2, C2, rq = res; fn = sys.argv[5]
json.dump({"N": N, "p": p, "q": q, "rankQ": rq, "E": E2, "C": C2}, open(fn, 'w')); print("rankQ", rq, "saved", fn)
