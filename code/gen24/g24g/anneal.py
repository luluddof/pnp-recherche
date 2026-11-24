# [E-search] simulated annealing over ternary E_l, C_l in {-1,0,1}^{N x N} minimising rank(Swap - sum E_l (x) C_l)
# (float SVD rank; any hit is re-verified exactly with flint by verify step at the end).
import numpy as np, sys, time, json, flint
N, q, steps, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
T0 = float(sys.argv[5]) if len(sys.argv) > 5 else 0.6
rng = np.random.default_rng(seed)
S = np.zeros((N*N, N*N))
for x in range(N):
    for y in range(N): S[x*N+y, y*N+x] = 1
def rk(E, C):
    K = sum(np.kron(e, c) for e, c in zip(E, C))
    s = np.linalg.svd(S - K, compute_uv=False)
    return int((s > 1e-8).sum()), float(np.log(s[s > 1e-8]).sum())
E = [rng.choice([-1.,0.,0.,0.,1.], size=(N,N)) for _ in range(q)]; C = [rng.choice([-1.,0.,0.,0.,1.], size=(N,N)) for _ in range(q)]
cur = rk(E, C); best = cur; bestEC = (E, C); t0 = time.time()
for it in range(steps):
    T = T0 * (1 - it/steps) + 1e-3
    l = rng.integers(q); side = rng.integers(2); i, j = rng.integers(N, size=2)
    M = (E if side == 0 else C)[l]; old = M[i, j]
    M[i, j] = rng.choice([v for v in (-1., 0., 1.) if v != old])
    new = rk(E, C)
    d = (new[0] - cur[0]) + 0.01*(new[1] - cur[1])
    if d <= 0 or rng.random() < np.exp(-d/T): cur = new
    else: M[i, j] = old
    if cur[0] < best[0]:
        best = cur; bestEC = ([e.copy() for e in E], [c.copy() for c in C])
        print(f"it {it} rank {best[0]} ({time.time()-t0:.0f}s)", flush=True)
E, C = bestEC
K = sum(np.kron(e, c) for e, c in zip(E, C)).astype(int); R = np.einsum('kxu,kyv->xuyv', np.array(E), np.array(C)).reshape(N*N, N*N).astype(int)
fz = lambda A: flint.fmpz_mat([[int(v) for v in r] for r in A])
rkQ = fz((S - K).astype(int)).rank(); osr = fz(R).rank()
print(f"FINAL N={N} q={q}: exact rank {rkQ}, exact osr {osr}, N(N-1)/2={N*(N-1)//2}")
json.dump({"N": N, "rank": rkQ, "osr": osr, "E": [e.astype(int).tolist() for e in E], "C": [c.astype(int).tolist() for c in C]}, open(f"anneal_{N}_{q}_s{seed}.json", "w"))
