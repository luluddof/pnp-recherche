# closed form of the (8,8)@8 Q8 array; quaternion units as (sign, letter)
import itertools, sys
from orbit import sign_solve, check_colouring, ev_check
sys.path.insert(0, '/home/user/pnp-recherche/code/gen22/g22b'); from shape import verify_shape
L = ['1', 'i', 'j', 'k']
tab = {('1', x): (1, x) for x in L}; tab.update({(x, '1'): (1, x) for x in L})
tab.update({('i', 'i'): (-1, '1'), ('j', 'j'): (-1, '1'), ('k', 'k'): (-1, '1'), ('i', 'j'): (1, 'k'), ('j', 'k'): (1, 'i'), ('k', 'i'): (1, 'j'),
            ('j', 'i'): (-1, 'k'), ('k', 'j'): (-1, 'i'), ('i', 'k'): (-1, 'j')})
U = [(s, x) for s in (1, -1) for x in L]
def mul(a, b): s, x = tab[(a[1], b[1])]; return (a[0] * b[0] * s, x)
def inv(a): return a if a[1] == '1' else (-a[0], a[1])
SB = {(-1, '1'), (1, 'i'), (1, 'j'), (1, 'k')}          # S-block: a transversal of Q8/{+-1}
tau = {'1': (1, '1'), 'i': (1, 'k'), 'j': (1, 'i'), 'k': (1, 'j')}
def colour(a, b):
    c = mul(inv(a), b); pal = 'S' if c in SB else 'T'
    return (pal, U.index(mul(a, tau[c[1]])) + (0 if pal == 'S' else 8))
col = [[colour(a, b) for b in U] for a in U]
print("colour violations:", check_colouring(col))
sg = sign_solve(col); print("signs solvable:", sg is not None)
nm = lambda u: ('' if u[0] > 0 else '-') + u[1]
print("sign s(a, a c) as table rows a, cols c:")
for a in U: print(' ', nm(a).rjust(3), ' '.join(('+' if sg[U.index(a)][U.index(mul(a, c))] > 0 else '-') for c in U))
cells = [[(col[x][y][0], col[x][y][1], sg[x][y]) for y in range(8)] for x in range(8)]
print(verify_shape(8, cells, quiet=True), ev_check(col, sg, dense=True))
from collections import Counter
print("class sizes:", Counter(Counter(col[x][y] for x in range(8) for y in range(8)).values()))
# invariant signs s(a, a c) = sigma(c)?  F_2 system in 8 unknowns
import numpy as np
rows = []
cls = {}
for x in range(8):
    for y in range(8): cls.setdefault(col[x][y], []).append((x, y))
for k, lst in cls.items():
    for (a, b), (a2, b2) in itertools.combinations(lst, 2):
        v = [0] * 8
        for (x, y) in ((a, b), (a2, b2), (a, b2), (a2, b)): v[U.index(mul(inv(U[x]), U[y]))] ^= 1
        rows.append(v + [1])
for sig in itertools.product((0, 1), repeat=8):
    if all((sum(r[i] * sig[i] for i in range(8)) % 2) == r[8] for r in rows): print("invariant sign sigma:", {nm(U[i]): '-' if sig[i] else '+' for i in range(8)}); break
else: print("no left-invariant sign sigma(c)")
