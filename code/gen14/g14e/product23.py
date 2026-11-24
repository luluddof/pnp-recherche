#!/usr/bin/env python3
"""2x3 grid, lam=1: product expansion from the one-bit identity [x=v][y=u] = d_xy d_uv - e_xu e_yv per rectangle;
grouped ranks for every threshold assignment."""
import itertools, numpy as np
from grid23_generic import T, rowflat, cutranks
d = np.eye(2); e = np.array([[0., 1.], [-1., 0.]])
# 12 bit axes: top j bit k (partner oth[j][k]) ; bottom j bit k.  axis order of T (dim-4 axes t1,b1,t2,b2,t3,b3), each = (bit0 low)
oth = {0: [1, 2], 1: [0, 2], 2: [0, 1]}
names = {}
for j in range(3):
    for k, j2 in enumerate(oth[j]):
        names[('t', j, j2)] = None; names[('b', j, j2)] = None
letters = iter('abcdefghijkl')
for key in names: names[key] = next(letters)
def out_sub():
    s = ''
    for j in range(3):
        for a in 'tb':
            s += names[(a, j, oth[j][1])] + names[(a, j, oth[j][0])]     # high bit first for reshape to dim 4
    return s
rects = [(0, 1), (0, 2), (1, 2)]
def term(choice):
    ops, subs = [], []
    for (j, j2), ch in zip(rects, choice):
        x, y, u, v = names[('t', j, j2)], names[('t', j2, j)], names[('b', j, j2)], names[('b', j2, j)]
        if ch == 'a': ops += [d, d]; subs += [x+y, u+v]
        else: ops += [-e, e]; subs += [x+u, y+v]
    return np.einsum(','.join(subs) + '->' + out_sub(), *ops).reshape((4,)*6)
terms = {c: term(c) for c in itertools.product('ab', repeat=3)}
S = sum(terms.values()); print("product identity sums to T:", np.abs(S - T).max() < 1e-12)
rr = lambda X: np.linalg.matrix_rank(rowflat(X), tol=1e-8)
for c, X in terms.items(): print("  term", ''.join(c), "row rank", rr(X), "col cut ranks", cutranks(X))
for k in range(0, 5):
    A = sum(X for c, X in terms.items() if c.count('b') < k) + 0*T; B = T - A
    print(f"  row group = terms with <{k} b's: row-rank(A)={rr(A)}, col cut ranks of B={cutranks(B)}; "
          f"sum-of-terms cost = {sum(4**c.count('b') if c.count('b')<k else 4**c.count('a') for c in terms)}")
x = np.load('lp_h3_x.npy'); print("LP (matching family) non-product solution: 6 terms, cost 24, grouped row-rank 10, col cut ranks [7,7,7]  (see grid23_h3.out)")
