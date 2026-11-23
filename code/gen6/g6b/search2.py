"""3- and 4-piece partitions of eps-matching representations (lambda0 = 1), ranked by the
multinomial grouped bound; then for the best, EXACT ranks of the type-class pieces for k = 2
blocks (lambda = 2), giving a rigorous lambda0 = 2 figure
   rho <= min_q log2 sum_n 2^{<q, c(G_n)>} / (4k)        (multinomial over the G_n themselves)."""
import itertools, math, random
import numpy as np
from box3 import costs, flat, rank
from skein import minimax, Tp
from grouped import rho_multinomial
import search1 as S1

random.seed(2)
cands = []
for key, a in S1.reps.items():
    sup = list(key)
    if len(sup) > 10:
        continue
    for npieces in (3, 4):
        n_all = npieces ** len(sup)
        labs = (itertools.product(range(npieces), repeat=len(sup)) if n_all <= 3000
                else (tuple(random.randrange(npieces) for _ in sup) for _ in range(3000)))
        seen = set()
        for lab in labs:
            if len(set(lab)) < npieces:
                continue
            groups = tuple(sorted(tuple(sup[j] for j in range(len(sup)) if lab[j] == g) for g in range(npieces)))
            if groups in seen:
                continue
            seen.add(groups)
            cs = [costs(S1.piece(a, list(g))) for g in groups]
            e, _ = minimax(cs)
            mu, q = rho_multinomial(cs)
            cands.append((mu, e / 4, key, groups, cs))
cands.sort(key=lambda r: r[0])
print("best by multinomial (3-4 pieces):")
for r in cands[:6]:
    print(f"  mult={r[0]:.4f} exp={r[1]:.4f} costs={[[round(x,2) for x in c] for c in r[4]]}")
cands.sort(key=lambda r: r[1])
print("best by single-block exp figure:")
for r in cands[:4]:
    print(f"  exp={r[1]:.4f} mult={r[0]:.4f} costs={[[round(x,2) for x in c] for c in r[4]]}")

# exact k=2 grouping for the top multinomial candidates
def kron8(A, B):
    # blocks: corner p gets the pair (w_p^A, w_p^B) -> index 2*a + b
    return np.einsum('abcdefgh,ijklmnop->aibjckdlemfngohp', A, B).reshape((4,) * 8)

cands.sort(key=lambda r: r[0])
print("exact k=2 type-class pieces:")
for r in cands[:3]:
    key, groups = r[2], r[3]
    a = S1.reps[key]
    P = [S1.piece(a, list(g)) for g in groups]
    s = len(P)
    types = {}
    for i, j in itertools.product(range(s), repeat=2):
        n = tuple(sorted((i, j)))
        types[n] = types.get(n, 0) + kron8(P[i], P[j])
    tot = sum(types.values())
    T2 = kron8(Tp, Tp)
    assert np.allclose(tot, T2)
    cs = [costs(G) for G in types.values()]
    mu, q = rho_multinomial(cs)
    e, _ = minimax(cs)
    print(f"  base mult={r[0]:.4f}: k=2 pieces={len(cs)} exp(k=2)={e/8:.4f} rigorous bound={mu*4/8:.4f}"
          f"  costs={[[round(x,2) for x in c] for c in cs]}")
