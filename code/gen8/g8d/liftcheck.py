"""d=1, generic tau: check  [G algebraically feasible]  <=>  [G balanced for gains (-1)^{|z|} e_k on edge z->z^k]
(i.e. every component of G lifts to the slab cover Lambda = {a in Z^l : sum a in {0,1}}, a ~ a+e_k for sum a = 0)."""
import random, sys
sys.path.insert(0, '.')
from lazy import Inst
def balanced(l, I, good):
    adj = {}
    for e in good:
        z, w, k = I.E[e]
        g = 1 if bin(z).count('1') % 2 == 0 else -1     # gain from z to w
        adj.setdefault(z, []).append((w, k, g)); adj.setdefault(w, []).append((z, k, -g))
    pot = {}
    for s in adj:
        if s in pot: continue
        pot[s] = (0,) * l; st = [s]
        while st:
            u = st.pop()
            for v, k, g in adj[u]:
                p = list(pot[u]); p[k] += g; p = tuple(p)
                if v in pot:
                    if pot[v] != p: return False
                else: pot[v] = p; st.append(v)
    return True
rng = random.Random(1); agree = tot = nfeas = 0
for l in (3, 4, 5):
    I = Inst(l, 1); ne = len(I.E)
    for trial in range(3000):
        p = rng.choice([0.3, 0.45, 0.55, 0.65])
        good = [e for e in range(ne) if rng.random() < p]
        a, b = I.feasible(good), balanced(l, I, good)
        tot += 1; agree += (a == b); nfeas += a
print(f"trials {tot}, agree {agree}, feasible {nfeas}")
