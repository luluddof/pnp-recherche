from orbit import *; from groups import catalogue
mul, inv, el = catalogue(8)['Q8']
# names: dicyclic(2): (i,s): a=(1,0) ~ 'i', x=(0,1) ~ 'j', a^2=-1
def nm(g):
    i, s = el[g]
    base = {0: '1', 1: 'i', 2: '-1', 3: '-i'}[i] if s == 0 else {0: 'j', 1: 'k', 2: '-j', 3: '-k'}[i]
    return base
# check (i,0)*(0,1) = (1,1) = 'k'? convention: (i,s): a^i x^s
names = [nm(g) for g in range(8)]
idx = {n: k for k, n in enumerate(names)}
assert mul[idx['i']][idx['j']] == idx['k'] and mul[idx['j']][idx['k']] == idx['i'] and mul[idx['k']][idx['i']] == idx['j']
E = Enc(mul, inv, [('S', frozenset([0])), ('T', frozenset([0]))], caps=[4, 4])
S = Cadical153(bootstrap_with=E.cl); sols = []
while S.solve():
    t, bl = E.decode(S.get_model()); sols.append((t, bl)); S.add_clause(E.block_clause(t, bl))
print("solutions with t_1 = 1:", len(sols))
from collections import Counter
sgn_ok = 0
for t, bl in sols:
    col = E.array(t, bl); sg = sign_solve(col); sgn_ok += sg is not None
print("sign-solvable:", sgn_ok)
Sblocks = Counter(tuple(sorted(names[c] for c in range(8) if bl[c] == 0)) for t, bl in sols)
print(Sblocks)
# pick solutions where S-block = {1,i,j,k} or similar; print r(c) = t_c^{-1} and base-class cells
for t, bl in sols[:400]:
    Sb = sorted(names[c] for c in range(8) if bl[c] == 0)
    if set(Sb) in ({'1', 'i', 'j', 'k'}, {'-1', 'i', 'j', 'k'}):
        print('S-block', Sb, ' t:', {names[c]: names[t[c]] for c in range(8)})
