# Exact test of Lemma 9 on a module with a DEGENERATE lowest-weight space on which zero modes act as a Clifford algebra:
# Ramond sector of n = 8 fermions, states a_x = b_x = psi^i psi^j (i<j; weight 1, so(8) currents, N = 28, H = 2),
# labels = the 16-dimensional ground space (E0 = 0), budget B = 2E0 + H = 2 (weights 0,1,2; dims 16, 128, 576).
# Rescaled modes {p^i_m, p^j_n} = 2 delta^{ij} delta_{m+n};  a[r] = sum_{m+n=-r} p^i_m p^j_n raises the weight by r.
import itertools, flint, sys
from ramond8 import G          # 8 real symmetric 16x16 gamma matrices (verified there)
n = 8
def lvl(mask): return bin(mask & 0xFF).count("1") + 2*bin(mask >> 8).count("1")
def app(i, m, vec):            # apply p^i_m to a dict {(mask, s): coef}
    out = {}
    if m == 0:
        for (mask, s), c in vec.items():
            sg = -1 if bin(mask).count("1") & 1 else 1
            for t in range(16):
                g = G[i][t][s]
                if g: out[(mask, t)] = out.get((mask, t), 0) + sg*g*c
        return out
    k = abs(m)
    if k > 2: return out if m > 0 else None      # annihilator kills level <= 2 ; creators beyond truncation never needed
    bit = 1 << (i + 8*(k-1))
    for (mask, s), c in vec.items():
        sg = -1 if bin(mask & (bit-1)).count("1") & 1 else 1
        if m < 0:
            if not mask & bit: out[(mask | bit, s)] = out.get((mask | bit, s), 0) + sg*c
        else:
            if mask & bit: out[(mask ^ bit, s)] = out.get((mask ^ bit, s), 0) + 2*sg*c
    return out
def J(i, j, r, vec):           # a[r] = sum_{m+nn=-r} p^i_m p^j_nn, annihilators applied first (i != j: they anticommute)
    out = {}
    for m in range(-2, 3):
        nn = -r - m
        if abs(nn) > 2: continue
        if m > 0: w = app(j, nn, app(i, m, vec)); sg = -1
        else: w = app(i, m, app(j, nn, vec)); sg = 1
        for k, c in w.items():
            if lvl(k[0]) <= 2: out[k] = out.get(k, 0) + sg*c
    return {k: c for k, c in out.items() if c}
pairs = list(itertools.combinations(range(n), 2)); N = len(pairs); Np = 16
ground = [(0, s) for s in range(16)]
# U[v][r][y] = b_v[r] |y>  (dict),  and the same vectors serve for <u| a_x[-r] via the adjoint: compute D by applying a[-r] to U
T = {s: {} for s in (-1, 0, 1)}
coef = {-1: {2: 1, 1: -2, 0: 1}, 0: {1: 1, 0: -2}, 1: {0: 1}}     # A_s = sum_r coef[s][r] * a[-r] b[r] on the ground space
for vi, (k, l) in enumerate(pairs):
    for yi, y in enumerate(ground):
        for r in (0, 1, 2):
            w = J(k, l, r, {y: 1})
            if not w: continue
            for xi, (i, j) in enumerate(pairs):
                res = J(i, j, -r, w)
                for (mask, u), c in res.items():
                    assert mask == 0
                    for s in coef:
                        if r in coef[s]:
                            key = (xi, yi, u, vi); T[s][key] = T[s].get(key, 0) + coef[s][r]*c
tot = {}
for s in T:
    for k, c in T[s].items(): tot[k] = tot.get(k, 0) + c
tot = {k: c for k, c in tot.items() if c}
ds = set(tot.values()); ok = all(k[0] == k[3] and k[1] == k[2] for k in tot) and len(tot) == N*Np and len(ds) == 1
print(f"Ramond n=8, so(8) currents (N={N}, H=2), labels = ground space (N'={Np}): sum_s A_s == d [x=v][y=u] exactly: {ok}, d = {ds}")
def rk(T, mode, xs=None):
    xs = list(range(N)) if xs is None else xs; pos = {x: a for a, x in enumerate(xs)}; n1 = len(xs)
    M = [[0]*(n1*Np) for _ in range(n1*Np)]
    for (x, y, u, v), c in T.items():
        if x in pos and v in pos:
            if mode == "r": M[pos[x]*Np + y][u*n1 + pos[v]] = c
            else: M[pos[x]*Np + u][y*n1 + pos[v]] = c
    return flint.fmpz_mat(M).rank()
def add(*Ts):
    o = {}
    for t in Ts:
        for k, c in t.items(): o[k] = o.get(k, 0) + c
    return {k: c for k, c in o.items() if c}
for name, xs in (("all 28 currents (28 x 16 rectangular crossing, N N' = 448)", None), ("first 16 currents (square N = 16, N^2 = 256)", list(range(16)))):
    print(name)
    for s in (-1, 0, 1): print(f"   piece s={s:2d}: rr = {rk(T[s], 'r', xs)}, rc = {rk(T[s], 'c', xs)}")
    print(f"   thresholds: (rr(A_-1), rc(A_0+A_1)) = ({rk(T[-1], 'r', xs)}, {rk(add(T[0], T[1]), 'c', xs)});  (rr(A_-1+A_0), rc(A_1)) = ({rk(add(T[-1], T[0]), 'r', xs)}, {rk(T[1], 'c', xs)})")
print("window dimensions: weight 0: 16, weights <=1: 144, weights <=2: 720")
