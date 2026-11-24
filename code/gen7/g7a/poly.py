"""Test algebraic feasibility of 'polychromatic' bad sets: edge (z,z^k) (z_k=0) bad iff
a=#ones of z before k == alpha mod p1 and b=#ones after k == beta mod p2 (AKS/Offner colouring classes),
plus layered classes (|z| == r mod m)."""
import sys
from core import edges, feasible
def test(l, d, badfn):
    E = edges(l); good = [e for e in E if not badfn(*e)]
    bad = len(E) - len(good)
    per = [sum(1 for e in E if e[2] == k and badfn(*e)) for k in range(l)]
    return feasible(l, d, good), bad, max(per)
for (l, d) in [(4,1),(5,1),(6,1),(4,2),(5,2),(6,2),(7,2),(5,3),(6,3),(7,3)]:
    m = d + 1
    res = []
    for p1 in range(1, m + 1):
        p2 = m + 1 - p1
        for al in range(p1):
            for be in range(p2):
                f = lambda z, w, k, p1=p1, p2=p2, al=al, be=be: (bin(z & ((1 << k) - 1)).count('1') % p1 == al and bin(z >> (k + 1)).count('1') % p2 == be)
                ok, bad, mx = test(l, d, f)
                if ok: res.append((p1, p2, al, be, bad, mx))
    lay = []
    for r in range(m):
        f = lambda z, w, k, r=r: bin(z).count('1') % m == r
        ok, bad, mx = test(l, d, f)
        lay.append(ok)
    print(f"l={l} d={d}: feasible polychromatic classes {res}  layered feasible {lay}", flush=True)
