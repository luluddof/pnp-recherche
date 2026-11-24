# "Matching ansatz" (pattern of all exact points of G19b): C = sum_l G_l (x) G'_l with G_l = sum_{(x,y) in M_l} s E_{x,y},
# G'_l = sum_{(x,y) in M_l} s E_{y,x}, M_l a partial permutation (distinct x, distinct y), the M_l disjoint, s = +-1.
# So C = 1 on the covered diagonal entries (x,y,y,x) and R = X - C must have rank <= p.  Simulated annealing on
# energy = rank(R) + w * (tail of singular values beyond p).  Found solutions are verified EXACTLY by verify.py.
import numpy as np, sys, time, json
from verify import verify
N, p, q = map(int, sys.argv[1].split(',')); tmax = float(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
rng = np.random.default_rng(seed)
X = np.zeros((N, N, N, N))
for x in range(N):
    for y in range(N): X[x, y, y, x] = 1
def build(assign, sign):
    G = np.zeros((q, N, N)); H = np.zeros((q, N, N))
    for x in range(N):
        for y in range(N):
            l = assign[x, y]
            if l >= 0: G[l, x, y] = sign[x, y]; H[l, y, x] = sign[x, y]
    return G, H
def energy(assign, sign):
    G, H = build(assign, sign)
    C = np.einsum('lxu,lyv->xyuv', G, H)
    s = np.linalg.svd((X - C).reshape(N*N, N*N), compute_uv=False)
    rk = int((s > 1e-9).sum())
    return max(rk - p, 0) + 0.02 * float((s[p:]**2).sum()), rk
def valid_move(assign, x, y, l):
    if l < 0: return True
    return not ((assign[x, :] == l).any() or (assign[:, y] == l).any())
best_overall = None; t0 = time.time(); restart = 0
while time.time() - t0 < tmax:
    restart += 1
    assign = -np.ones((N, N), int); sign = rng.choice([-1, 1], size=(N, N))
    # random greedy init: cover about half of the pairs
    for x, y in rng.permutation([(a, b) for a in range(N) for b in range(N)]):
        if rng.random() < 0.6:
            for l in rng.permutation(q):
                if valid_move(assign, x, y, l): assign[x, y] = l; break
    E, rk = energy(assign, sign); T = 1.0; it = 0; bestE = E
    while time.time() - t0 < tmax and it < 40000:
        it += 1
        x, y = rng.integers(N), rng.integers(N); old = (assign[x, y], sign[x, y])
        r = rng.random()
        if r < 0.3: sign[x, y] *= -1
        else:
            l = rng.integers(-1, q)
            if l == assign[x, y] or not valid_move(assign, x, y, l): continue
            assign[x, y] = l
        E2, rk2 = energy(assign, sign)
        if E2 <= E or rng.random() < np.exp((E - E2) / T):
            E, rk = E2, rk2
        else:
            assign[x, y], sign[x, y] = old
        T = max(0.02, T * 0.9997)
        if rk <= p: break
        if it % 5000 == 0: print(f"  restart {restart} it {it} E {E:.3f} rank {rk} T {T:.3f}", flush=True)
    print(f"restart {restart}: final rank {rk} (target {p}), E {E:.3f}, its {it}, time {time.time()-t0:.0f}s", flush=True)
    if rk <= p:
        G, H = build(assign, sign)
        pp, qq = verify(N, G.astype(int).tolist(), H.astype(int).tolist(), f"matching ansatz N={N}")
        json.dump({"N": N, "p": pp, "q": qq, "G": [[[str(int(e)) for e in r] for r in l] for l in G], "Gp": [[[str(int(e)) for e in r] for r in l] for l in H]}, open(f"exact_m_{N}_{pp}_{qq}.json", "w"))
        print("FOUND exact", (pp, qq)); break
