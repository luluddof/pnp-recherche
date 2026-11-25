# Exact (over Z/Q) feasibility check of the degree-1 ansatz, D = h, no sketch, for small h.
import sys, time
sys.dont_write_bytecode = True
import flint
from math import comb
from appb_exact import partitions, apply_mono

h = int(sys.argv[1]); tests = [tuple(int(a) for a in t.split(',')) for t in sys.argv[2:]]   # (W, W')
bas = {}; idx = {}
def B(E):
    if E not in bas:
        bas[E] = sorted(partitions(E)) if E >= 0 else []; idx[E] = {s: i for i, s in enumerate(bas[E])}
    return bas[E]
def mat(mu, r, E):
    src = B(E); tgt = B(E + r); M = flint.fmpz_mat(max(len(tgt), 1), max(len(src), 1))
    for j, s in enumerate(src):
        for s2, c in apply_mono(mu, r, s): M[idx[E + r][s2], j] += c
    return M
Vh = B(h); n = len(Vh)
mus = [m for w in range(h + 1) for m in sorted(partitions(w))]
cols = {}
for mu in mus:
    hp = sum(mu); H = h + hp; fp = {}; fm = {}
    for w in range(0, 2 * h + 1):
        if not B(w): continue
        Ad = mat(mu, h - w, w); Au = mat(mu, w - h, h)
        fpv = []; fmv = []
        for nu in Vh:
            X = Ad * mat(nu, w - h, h)        # t^mu[h-w] t^nu[w-h] on V_h
            Y = mat(nu, h - w, w) * Au        # t^nu[h-w] t^mu[w-h] on V_h
            fpv += [int(X[i, j]) for i in range(n) for j in range(n)]
            fmv += [int(Y[i, j]) for i in range(n) for j in range(n)]
        fp[w - 2 * h] = fpv; fm[-w] = fmv
    L = n ** 3
    def Fe(e, f):
        acc = [0] * L
        for j in range(H + 1):
            if (e - j) in f:
                cj = (-1) ** j * comb(H, j); v = f[e - j]
                acc = [a + cj * b for a, b in zip(acc, v)]
        return acc
    for e in range(-2 * h, h + hp + 1):
        cols[(mu, e)] = Fe(e, fp) if e <= 0 else Fe(e, fm)
rhs = []
for k0 in range(n):
    v = [0] * n ** 3
    for s in range(n): v[k0 * n * n + s * n + s] = 1
    rhs.append(v)
def rank(cl):
    M = flint.fmpz_mat(len(cl[0]), len(cl), [cl[j][i] for i in range(len(cl[0])) for j in range(len(cl))])
    return M.rank()
print(f"h={h}: {len(cols)} columns, exact integer arithmetic; rank of all = {rank(list(cols.values()))}", flush=True)
for W, Wp in tests:
    cl = [v for (mu, e), v in cols.items() if e <= Wp - 2 * h or e >= sum(mu) + h - W]
    ra, rb = rank(cl), rank(cl + rhs)
    print(f"  W={W} W'={Wp} k={4*h-1-W-Wp}: rank A={ra}, rank [A|Id(x)l]={rb} -> feasible over Q: {ra == rb}", flush=True)
