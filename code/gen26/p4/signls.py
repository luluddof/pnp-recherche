# local search over signs for a fixed support (objective: float rank / sum of small singular values), exact check at the end
import sys, json, numpy as np, mc, time
N, p, q = map(int, sys.argv[1:4]); solf = sys.argv[4]; out = sys.argv[5]; TL = float(sys.argv[6])
ms = set()
for line in open(solf):
    if line.startswith('v'): ms.update(int(t) for t in line.split()[1:] if int(t) > 0)
idx = lambda l, x, u, side: 1 + side * q * N * N + l * N * N + x * N + u
Es = np.array([[[int(idx(l, x, u, 0) in ms) for u in range(N)] for x in range(N)] for l in range(q)])
Cs = np.array([[[int(idx(l, y, v, 1) in ms) for v in range(N)] for y in range(N)] for l in range(q)])
pos = [(0, l, x, u) for l in range(q) for x in range(N) for u in range(N) if Es[l, x, u]] + \
      [(1, l, x, u) for l in range(q) for x in range(N) for u in range(N) if Cs[l, x, u]]
S = mc.swap(N).astype(float); rng = np.random.default_rng(0); t0 = time.time()
def build(sg):
    E = Es.astype(float).copy(); C = Cs.astype(float).copy()
    for s, (side, l, x, u) in zip(sg, pos):
        (E if side == 0 else C)[l, x, u] = s
    return E, C
def f(sg):
    E, C = build(sg); K = sum(np.kron(E[l], C[l]) for l in range(q))
    sv = np.linalg.svd(S - K, compute_uv=False); return sv[p:].sum()
while time.time() - t0 < TL:
    sg = rng.choice([-1., 1.], len(pos)); cur = f(sg)
    for it in range(3000):
        i = rng.integers(len(pos)); sg[i] *= -1; v = f(sg)
        if v <= cur + 1e-12: cur = v
        else: sg[i] *= -1
        if cur < 1e-9: break
    if cur < 1e-9:
        E, C = build(sg); E = E.astype(int).tolist(); C = C.astype(int).tolist()
        rq = mc.rankQ(N, E, C); print("found rankQ", rq, f"[{time.time()-t0:.1f}s]")
        if rq <= p: json.dump({"N": N, "p": p, "q": q, "rankQ": rq, "E": E, "C": C}, open(out, 'w')); print("saved", out); break
else: print("no signing found in", TL, "s")
