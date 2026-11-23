"""Search at lambda0 = 1 over pieces built from eps-matching representations of T'.
1. vertex representations T' = sum a_M eps_M (weighted L1 minimisation, random weights)
2. for each: every 2-partition of the support -> pieces (A, B); figures:
   exp  = min_q max_s <q,c(s)>/4            (single block; exponentially many terms)
   mult = multinomial upper bound, grouped  (poly pieces)
   grp2 = exact eigen-count grouped figure (2 pieces)
3. 3-partitions sampled at random -> exp and mult figures."""
import itertools, math, random, sys
import numpy as np
from scipy.optimize import linprog
from box3 import all_matchings, CORNERS, costs, flat
from skein import eps_matching, cross, Tp, minimax
from grouped import rho_two_piece, rho_multinomial

Ms = list(all_matchings(CORNERS))
EM = np.stack([eps_matching(M).ravel() for M in Ms], 1)      # 256 x 105
t = Tp.ravel()

def l1_rep(w):
    n = len(Ms)
    # a = a+ - a- ; min w.(a+ + a-) ; EM (a+ - a-) = t
    res = linprog(np.concatenate([w, w]), A_eq=np.hstack([EM, -EM]), b_eq=t,
                  bounds=[(0, None)] * (2 * n), method="highs")
    a = res.x[:n] - res.x[n:]
    a[abs(a) < 1e-9] = 0
    return a

def piece(a, idx):
    return (EM[:, idx] @ a[idx]).reshape((2,) * 8)

random.seed(1); np.random.seed(1)
reps = {}
base_w = np.array([2.0 ** sum(cross(M)) for M in Ms])
for trial in range(60):
    w = base_w * np.exp(np.random.randn(len(Ms)) * (0 if trial == 0 else 0.7))
    a = l1_rep(w)
    key = tuple(np.nonzero(a)[0])
    reps[key] = a
print("distinct vertex representations:", len(reps), " support sizes:", sorted(len(k) for k in reps)[:10])

results = []
for key, a in reps.items():
    sup = list(key)
    if len(sup) > 16:
        continue
    single = [tuple(round(x, 3) for x in costs(piece(a, [m]))) for m in sup]
    exp1, _ = minimax(single)
    for mask in range(1, 2 ** (len(sup) - 1)):
        I1 = [sup[j] for j in range(len(sup)) if mask >> j & 1]
        I2 = [m for m in sup if m not in I1]
        A, B = piece(a, I1), piece(a, I2)
        cA, cB = costs(A), costs(B)
        e, _ = minimax([cA, cB])
        if e / 4 > 0.75:
            continue
        mu, _ = rho_multinomial([cA, cB])
        results.append((e / 4, mu, key, mask, cA, cB))
results.sort(key=lambda r: (r[1], r[0]))
print("best 2-partitions by multinomial grouped figure:")
for r in results[:8]:
    print(f"  exp={r[0]:.4f} mult={r[1]:.4f} cA={[round(x,2) for x in r[4]]} cB={[round(x,2) for x in r[5]]}")
# exact eigen-count grouped figure for top candidates
print("exact 2-piece grouped figure (eigen count) for top 30 by multinomial:")
best = None
for r in results[:30]:
    key, mask = r[2], r[3]
    a = reps[key]; sup = list(key)
    I1 = [sup[j] for j in range(len(sup)) if mask >> j & 1]
    A = piece(a, I1)
    rg, q, ok, F = rho_two_piece(A, Tp, 1.0)
    if best is None or rg < best[0]:
        best = (rg, q, ok, r)
    print(f"  grp={rg:.4f} diag={ok} exp={r[0]:.4f} mult={r[1]:.4f}")
print("BEST exact grouped 2-piece:", best[0], "q=", best[1], "diagonalizable:", best[2])
sys.stdout.flush()
