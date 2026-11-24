# G10c: is the encoding really O(m log p) clauses of width 3 (no hidden log m from the cardinality constraints)?
import random, math, sys
sys.path.insert(0, '.')
from encode import build, CNF
rnd = random.Random(5)
print("counter alone: t, gates/t, clauses/t")
for t in (10, 100, 1000, 10000, 100000):
    F = CNF(); bits = [F.new() for _ in range(t)]; o = F.counter(bits); F.equals_const(o, t // 3)
    print(f"  {t:8d} {F.gates / t:.3f} {len(F.cl) / t:.3f}  output bits {len(o)}")
print("full Psi (Steiner / cycle / path), random 3-regular-ish graphs, m = 1.5 n: kind, m, p, b, clauses/(b*m), clauses/m, vars/m")
for kind in ('ST', 'CYC', 'PATH'):
    for n in (200, 2000, 20000):
        E = set()
        while len(E) < 3 * n // 2:
            u, v = rnd.randrange(n), rnd.randrange(n)
            if u != v: E.add((min(u, v), max(u, v)))
        E = sorted(E); m = len(E)
        for p in (2, 5, 17, 257):
            w = [rnd.randint(1, p - 1) if p > 2 else 1 for _ in E]
            F, M0, b = build(n, E, p, w, kind, n // 2, K=[0, 1, 2])
            assert max(len(c) for c in F.cl) <= 3
            print(f"  {kind:4s} m={m:7d} p={p:3d} b={b} cl/(b m)={len(F.cl) / (b * m):6.2f} cl/m={len(F.cl) / m:7.2f} vars/m={F.nv / m:6.2f}  skeleton cl/m={M0 / m:5.2f}")
