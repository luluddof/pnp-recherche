#!/usr/bin/env python3
"""G18S r1: independent recomputation of Cor A'1 on the uniform product-label instance.
Instance U(u,r): labels [r]^u, ONE wire per ordered pair (row label x, column label y): M = r^(2u).
Direction i: group = symbol x_i.  Everything is computed by explicit enumeration of wires (no formula),
then compared with the closed forms."""
import itertools, math
from fractions import Fraction as Fr

def aprime_vertex_cut(u, r, inC1, inC2):
    """vertex(label)-cut version (the notes' A'). returns (w, [c_i]) by explicit enumeration."""
    labels = list(itertools.product(range(r), repeat=u))
    cut = []  # (label of C1-end, label of C2-end)
    for x in labels:
        for y in labels:
            if inC1(x) and inC2(y): cut.append((x, y))
            elif inC2(x) and inC1(y): cut.append((y, x))
    w = len(cut); cs = []
    for i in range(u):
        e1 = [0] * r; e2 = [0] * r
        for (p, q) in cut: e1[p[i]] += 1; e2[q[i]] += 1
        cs.append(sum(min(a, b) for a, b in zip(e1, e2)))
    return w, cs

def aprime_oriented(u, r, delta):
    """generalised (oriented wire subset) version A'': delta(x,y) in {+1: row end is the 'row index',
    -1: column end is, 0: wire not selected}."""
    labels = list(itertools.product(range(r), repeat=u))
    sel = []
    for x in labels:
        for y in labels:
            d = delta(x, y)
            if d > 0: sel.append((x, y))
            elif d < 0: sel.append((y, x))
    w = len(sel); cs = []
    for i in range(u):
        e1 = [0] * r; e2 = [0] * r
        for (p, q) in sel: e1[p[i]] += 1; e2[q[i]] += 1
        cs.append(sum(min(a, b) for a, b in zip(e1, e2)))
    return w, cs

print("u  r | majority vertex cut: w/M, (w-max c)/M | closed form | notes' value | oriented A'': (w-max c)/M | C(2u,u)/4^u")
for r in (2, 4):
    for u in range(1, 13):
        if r ** (2 * u) > 2 ** 20 * 20: continue
        half = lambda s: 1 if s < r // 2 else 0
        sc = lambda x: sum(half(s) for s in x)
        if u % 2:
            C1 = lambda x: sc(x) * 2 > u; C2 = lambda x: sc(x) * 2 < u
        else:   # even u: A = {s >= u/2}, B = {s < u/2}  (the notes' best threshold)
            C1 = lambda x: sc(x) * 2 >= u; C2 = lambda x: sc(x) * 2 < u
        M = r ** (2 * u)
        w, cs = aprime_vertex_cut(u, r, C1, C2)
        uu = u if u % 2 else u + 1
        closed = Fr(math.comb(uu - 1, (uu - 1) // 2), 2 ** uu)
        w2, cs2 = aprime_oriented(u, r, lambda x, y: (sc(x) > sc(y)) - (sc(x) < sc(y)))
        print(f"{u:2d} {r:2d} | {Fr(w, M)} {Fr(w - max(cs), M)} = {float(Fr(w - max(cs), M)):.5f} | {float(closed):.5f} | "
              f" | {Fr(w2 - max(cs2), M)} = {float(Fr(w2 - max(cs2), M)):.5f} | {math.comb(2*u, u)/4**u:.5f}")
print("closed forms for larger u (majority: C(u'-1,(u'-1)/2)/2^u', u' = u or u+1 odd; oriented: C(2u,u)/4^u):")
for u in list(range(2, 13)) + [41, 401]:
    uu = u if u % 2 else u + 1
    a = math.comb(uu - 1, (uu - 1) // 2) / 2 ** uu; b = math.comb(2 * u, u) / 4 ** u
    print(f"  u={u:3d}  majority {a:.5f}  (*sqrt u = {a*math.sqrt(u):.4f})   oriented {b:.5f} (*sqrt u = {b*math.sqrt(u):.4f})   Thm A {2/4**u:.2e}")
print("limits: 1/sqrt(2 pi) =", 1 / math.sqrt(2 * math.pi), "  1/sqrt(pi) =", 1 / math.sqrt(math.pi))
