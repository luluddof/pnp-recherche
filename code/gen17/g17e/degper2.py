from degper import *
import sympy as sp
# reproduce the (3,1) layered merged case with deg 0, symbolically over Q(tau)
rng = random.Random(3)
l, d = 3, 1; n = 8
# replay rng consumption of degper.main: runs (3,1) layered full, select, merged
tau = [rng.randrange(2, P - 1) for _ in range(l)]            # full
tau = [rng.randrange(2, P - 1) for _ in range(l)]; U = rng.sample(range(n), n // 2)   # select
tau = [rng.randrange(2, P - 1) for _ in range(l)]
a = list(range(n)); rng.shuffle(a); b = list(range(n)); rng.shuffle(b)
rows = [a[2 * i:2 * i + 2] for i in range(4)]; cols = [b[2 * i:2 * i + 2] for i in range(4)]
print("rows", rows, "cols", cols, "blocks", blocks_layered(l, d))
t = sp.symbols("t"); b_ = sp.symbols("b1:4")   # b_j = 1/tau_j
bl = blocks_layered(l, d); g = []
for z in range(n):
    w = bin(z).count("1"); s = [1 - 2 * (z >> j & 1) for j in range(l)]
    B = next(B for B in bl if w in B)
    Uv = sum(s[j] * (1 - t * b_[j]) for j in range(l)); x = sp.symbols("x")
    p = sp.interpolate([(l - 2 * ww, (-1) ** ww) for ww in B], x)
    g.append(sp.expand((-1) ** w + p.subs(x, Uv) - p.subs(x, sum(s))))
gam = [sp.expand(sum(g[z] * (-1) ** bin(z & w).count("1") for z in range(n)) / n) for w in range(n)]
A = sp.Matrix(4, 4, lambda i, j: sp.expand(sum(gam[x ^ y] for x in rows[i] for y in cols[j])))
sp.pprint(A); print("per =", sp.factor(sp.expand(A.per())))
print("gamma_w:", {w: gam[w] for w in range(n)})
