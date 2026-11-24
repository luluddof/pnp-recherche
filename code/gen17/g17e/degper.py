"""G18d Task B: degree of P(t) = per(L J(t) R') for the layered / single-F families (commutative class), mod p.
Also the 2x2 cancellation examples (sympy)."""
import random, itertools, numpy as np, sympy as sp
P = 2147483629
def inv(a): return pow(a % P, P - 2, P)
def per_mod(A):
    n = len(A); A = np.array(A, dtype=np.int64) % P
    S = np.zeros((1, n), dtype=np.int64); sign = np.ones(1, dtype=np.int64)
    for i in range(n):                       # subsets of ROWS; row sums vector over columns
        S = np.concatenate([S, (S + A[i]) % P]); sign = np.concatenate([sign, (-sign) % P])
    prod = np.ones(len(S), dtype=np.int64)
    for j in range(n): prod = prod * S[:, j] % P
    tot = int((prod * sign % P).sum() % P)
    return tot * (1 if n % 2 == 0 else P - 1) % P
def interp_deg(xs, ys):
    """degree of the interpolating polynomial through (xs, ys) mod P (Newton divided differences)."""
    n = len(xs); c = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            c[i] = (c[i] - c[i - 1]) * inv(xs[i] - xs[i - j]) % P
    deg = max([i for i in range(n) if c[i] != 0], default=-1)
    return deg
def pB_coeffs(l, block):
    """p with p(l-2w) = (-1)^w, w in block; returns function evaluating p at a value mod P (Lagrange)."""
    pts = [(l - 2 * w) % P for w in block]; vals = [1 if w % 2 == 0 else P - 1 for w in block]
    def ev(s):
        tot = 0
        for i, (a, v) in enumerate(zip(pts, vals)):
            num = 1; den = 1
            for j, b in enumerate(pts):
                if j != i: num = num * ((s - b) % P) % P; den = den * ((a - b) % P) % P
            tot = (tot + v * num % P * inv(den)) % P
        return tot
    return ev
def blocks_layered(l, d):
    # central block of d+1 weights, then blocks of d+1 outward
    w0 = (l - d) // 2; cuts = sorted(set([0, l + 1] + [w0 + i * (d + 1) for i in range(-l, l) if 0 < w0 + i * (d + 1) <= l]))
    return [list(range(cuts[i], cuts[i + 1])) for i in range(len(cuts) - 1)]
def J_at(l, d, tau, t, family):
    """J(t) as 2^l x 2^l integer matrix mod P, J = 2^{-l} H diag(g_z(t)) H."""
    n = 1 << l
    if family == "layered": bl = blocks_layered(l, d)
    else: w0 = (l - d) // 2; bl = [list(range(w0, w0 + d + 1))]
    ps = [pB_coeffs(l, b) for b in bl]
    g = []
    for z in range(n):
        w = bin(z).count("1"); s = [1 - 2 * (z >> j & 1) for j in range(l)]
        if family == "layered": bi = next(i for i, b in enumerate(bl) if w in b)
        else: bi = 0
        U = sum(s[j] * (1 - t * inv(tau[j])) for j in range(l)) % P
        chi = 1 if w % 2 == 0 else P - 1
        g.append((chi + ps[bi](U) - ps[bi](sum(s) % P)) % P)
    i2 = inv(n); J = [[0] * n for _ in range(n)]
    gam = [sum((g[z] if bin(z & w).count("1") % 2 == 0 else P - g[z]) for z in range(n)) % P * i2 % P for w in range(n)]
    for x in range(n):
        for y in range(n): J[x][y] = gam[x ^ y]
    return J
def run(l, d, family, inp, rng):
    n = 1 << l; tau = [rng.randrange(2, P - 1) for _ in range(l)]
    if inp == "full": rows = [[x] for x in range(n)]; cols = [[y] for y in range(n)]
    elif inp == "select":
        U = rng.sample(range(n), n // 2); rows = [[x] for x in U]; cols = [[(n - 1) ^ x] for x in U]
    else:  # merged: n/2 rows, each = 2 cube row-vertices; columns = 2 column-vertices (random pairing)
        a = list(range(n)); rng.shuffle(a); b = list(range(n)); rng.shuffle(b)
        rows = [a[2 * i:2 * i + 2] for i in range(n // 2)]; cols = [b[2 * i:2 * i + 2] for i in range(n // 2)]
    m = len(rows); xs = list(range(1, m * d + 3)); ys = []
    def A_of(J): return [[sum(J[x][y] for x in R for y in C) % P for C in cols] for R in rows]
    J0 = J_at(l, d, tau, 0, family)
    assert all(J0[x][y] == (1 if y == (n - 1) ^ x else 0) for x in range(n) for y in range(n))
    for t in xs: ys.append(per_mod(A_of(J_at(l, d, tau, t, family))))
    dg = interp_deg(xs, ys); p0 = per_mod(A_of(J0))
    print(f"  l={l} d={d} {family:8s} input={inp:6s}: n={m}  n*d={m*d}  deg per A(t) = {dg}   per A(0)={p0 if p0 < P//2 else p0-P}", flush=True)
if __name__ == "__main__":
    s, T = sp.symbols("s T")
    J1 = sp.Matrix([[s, 1 - s], [1 - s, -s]]); print("2x2, d=1: per =", sp.expand(sp.Matrix.per(J1)), " (n d = 2)")
    r2 = sp.sqrt(2)
    J2 = sp.Matrix([[s * (r2 - s), 1 - s**2], [1 - s**2, s * (r2 + s)]]); print("2x2, d=2: per =", sp.simplify(sp.expand(sp.Matrix.per(J2))), " (n d = 4)")
    print("   J2(s=0) =", J2.subs(s, 0).tolist(), " J2(s=1) =", [[sp.simplify(e) for e in r] for r in J2.subs(s, 1).tolist()])
    rng = random.Random(3)
    for l, d in [(3, 1), (3, 2), (4, 1), (4, 2), (4, 3)]:
        for family in ["layered", "singleF"]:
            for inp in ["full", "select", "merged", "merged"]:
                if l == 4 and inp == "full" and d == 3: continue
                run(l, d, family, inp, rng)
