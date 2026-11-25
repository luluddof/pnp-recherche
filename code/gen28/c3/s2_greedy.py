# Greedy randomized search for max mu_{1/4}-agreement sets Z of degree<=d multilinear g with parity (site-pattern model).
# Consistency maintained exactly mod p by incremental echelon form on augmented rows [A_z | par_z].
import numpy as np, itertools, sys
from math import comb
PR = 2147483629
class Ech:
    def __init__(s, M): s.rows = {}; s.M = M   # pivot col -> row (normalized)
    def reduce(s, v):
        v = v.copy()
        for c in sorted(s.rows):
            if v[c]: v = (v - v[c] * s.rows[c]) % PR
        return v
    def insert(s, v):
        c = int(np.nonzero(v)[0][0]); inv = pow(int(v[c]), PR-2, PR); v = (v * inv) % PR
        for c2 in list(s.rows):
            r = s.rows[c2]
            if r[c]: s.rows[c2] = (r - r[c] * v) % PR
        s.rows[c] = v
def run(k, d, restarts, seed=0):
    rng = np.random.default_rng(seed)
    pats = np.array(list(itertools.product([0,1], repeat=k)))
    mu = np.prod(np.where(pats == 1, 0.25, 0.75), axis=1)
    mons = [U for r in range(d+1) for U in itertools.combinations(range(k), r)]
    M = len(mons)
    A = np.array([[int(np.prod(T[list(U)])) if U else 1 for U in mons] + [(-1)**int(T.sum()) % PR] for T in pats], dtype=np.int64)
    best = 0; bestZ = None
    for r in range(restarts):
        noise = rng.random(len(pats)) ** (0.0 if r == 0 else rng.uniform(0.2, 3))
        order = np.argsort(-(mu * noise))
        E = Ech(M+1); Z = []; tot = 0.0
        for i in order:
            v = E.reduce(A[i] % PR)
            nz = np.nonzero(v)[0]
            if len(nz) == 0: Z.append(i); tot += mu[i]; continue
            if nz[0] == M: continue      # would make system inconsistent (pivot in parity column)
            E.insert(v); Z.append(i); tot += mu[i]
        if tot > best + 1e-12: best = tot; bestZ = sorted(Z)
    return best, bestZ, pats, mu
if __name__ == '__main__':
    from s2_lp import sym_best, bound
    k = int(sys.argv[1]); R = int(sys.argv[2])
    for d in range(1, int(sys.argv[3]) if len(sys.argv) > 3 else k):
        b, Z, pats, mu = run(k, d, R, seed=d)
        # symmetry check of best Z: is it a union of weight classes?
        ws = {}
        for i in Z: ws[int(pats[i].sum())] = ws.get(int(pats[i].sum()), 0) + 1
        full = {m: comb(k, m) for m in range(k+1)}
        prof = ' '.join(f"{m}:{ws.get(m,0)}/{full[m]}" for m in range(k+1))
        print(f"k={k} d={d} greedy={b:.4f} sym={sym_best(k,d):.4f} bound={bound(k,d):.4f}  weight profile {prof}", flush=True)
