# Check the two lemmas behind Theorem S' at d=2:
# (L1) span{ pi (X^{(x)(k-u)} (x) Y) pi^-1 } = sum of conjugation-isotypic components mu with mu_1 >= k-u (dimension check)
# (L2) such operators map total-spin sector m' (shape (k-m',m')) only into sectors m with |m-m'| <= u.
import numpy as np, itertools
from math import comb, factorial
rng = np.random.default_rng(2)

def perm_op(k, perm, d=2):
    # operator permuting tensor factors of (C^d)^{(x)k}: e_{i_1..i_k} -> e_{i_{perm^-1}}
    n = d**k; M = np.zeros((n,n))
    for idx in itertools.product(range(d), repeat=k):
        new = [None]*k
        for s in range(k): new[perm[s]] = idx[s]
        a = int(np.ravel_multi_index(idx, [d]*k)); b = int(np.ravel_multi_index(new, [d]*k))
        M[b,a] = 1
    return M

def kron_list(ms):
    out = np.array([[1.0]])
    for m in ms: out = np.kron(out, m)
    return out

def spin_sector_projectors(k):
    # total spin J = k/2 - m; Casimir J^2 = sum of S_i.S_j ; use sum of swaps T = sum_{i<j} SWAP_ij, eigenvalue on shape (k-m,m):
    # content sum = C(k-m,2) - (m choose 2)... compute: contents row1: 0..k-m-1, row2: -1..m-2
    n = 2**k; T = np.zeros((n,n))
    for i in range(k):
        for j in range(i+1,k):
            perm = list(range(k)); perm[i], perm[j] = j, i
            T += perm_op(k, perm)
    w, V = np.linalg.eigh(T)
    projs = {}
    for m in range(k//2+1):
        c = sum(range(k-m)) + sum(range(-1, m-1))
        sel = np.abs(w - c) < 1e-6
        projs[m] = V[:, sel] @ V[:, sel].T
    return projs

def rand_gen_op(k, u):
    X = rng.normal(size=(2,2)); Y = rng.normal(size=(2**u, 2**u))
    op = np.kron(kron_list([X]*(k-u)), Y)
    perm = list(rng.permutation(k)); Pm = perm_op(k, perm)
    return Pm @ op @ Pm.T

for k in [4,5,6]:
    projs = spin_sector_projectors(k)
    for u in range(0, k//2+1):
        viol = 0.0
        for t in range(20):
            A = sum(rng.normal()*rand_gen_op(k,u) for _ in range(5))
            for m in projs:
                for mp in projs:
                    if abs(m-mp) > u:
                        viol = max(viol, np.abs(projs[m] @ A @ projs[mp]).max())
        # also allowed blocks are generically nonzero:
        A = rand_gen_op(k,u); allowed = {(m,mp): np.abs(projs[m]@A@projs[mp]).max() > 1e-8 for m in projs for mp in projs}
        print(f"k={k} u={u}: max forbidden-block entry = {viol:.2e}; nonzero blocks |dm|: {sorted(set(abs(a-b) for (a,b),v in allowed.items() if v))}")
# dimension check of L1 at k=4 (labels = End(C^2)^{(x)4}, 256-dim)
k = 4
def isotypic_dims(k, D):
    # sum_{mu_1 >= k-u} f^mu * dim S^mu(C^D) via hook formulas
    from sympy.utilities.iterables import partitions
    res = {}
    for p in partitions(k):
        mu = sorted([a for a,c in p.items() for _ in range(c)], reverse=True)
        if len(mu) > D: continue
        cols = [sum(1 for r in mu if r > c) for c in range(mu[0])]
        hooks = 1; contprod = 1
        for i,r in enumerate(mu):
            for c in range(r):
                hooks *= (r - c - 1) + (cols[c] - i - 1) + 1
                contprod *= (D + c - i)
        f = factorial(k)//hooks; dimW = contprod//hooks
        res[tuple(mu)] = f*dimW
    return res
dims = isotypic_dims(k, 4)
for u in range(0, k):
    pred = sum(v for mu,v in dims.items() if mu[0] >= k-u)
    vecs = np.array([rand_gen_op(k,u).ravel() for _ in range(400)])
    r = np.linalg.matrix_rank(vecs, tol=1e-7)
    print(f"L1 dim check k={k} u={u}: span dim (400 samples) = {r}, predicted sum_(mu1>=k-u) = {pred}")
