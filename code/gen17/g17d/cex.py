# counterexamples to g14a sec.3 when the nonsingular-minor hypothesis is dropped
import itertools
from fractions import Fraction
from math import comb, ceil
import numpy as np
from test1 import rank, binp
P = 101
# (a) q=2, l=6, d=0: A(t) = all-ones 64x64 (constant), alpha=beta=id, pi = complement, B_k=0, E_k=A.
l, d = 6, 0; m = (l - d - 1) // 2; n = 2 ** l
A = np.ones((n, n), dtype=np.int64)
tot = l * rank(A, P)
print("(a) sum rank E_k =", tot, " claimed bound n p^2 =", float(n * binp(l, m) ** 2), "-> violated:", tot < n * binp(l, m) ** 2)
# (b) every nonzero has property (iii): q=2, l=1, n=5, all rows label 0, all cols label 1, A=all-ones, d=0.
l, d, n = 1, 0, 5; m = 0
print("(b) rank E_1 = 1, claimed bound n p^2 =", float(n * binp(l, m) ** 2), "-> violated:", 1 < n * binp(l, m) ** 2)
# (c) injective labels + property (iii) for every nonzero: q=4, rows {0,1}^l, cols {2,3}^l, all-ones, l=6
l, d = 6, 0; m = (l - d - 1) // 2; n = 2 ** l
print("(c) q=4 l=6: sum rank E_k =", l, " claimed", float(n * binp(l, m) ** 2))
# (d) sparse 0/1 A with per(A) != 0: n/2 blocks [[1,1],[1,1]]; rows/cols of block b labelled so that all 4 entries satisfy (iii):
#     rows: words w_b, w_b' ; cols: words differing from both in every coordinate (q=4). A(t)=A, E_k=A, sum rank = l*n/2
#     vs n p^2 (p ~ 1/2 at d=0): l n/2 >= n/4, no violation -- rank deficiency of a Delta-sparse A costs only a constant.
print("(d) 2x2 all-ones blocks: sum rank E_k = l*n/2 >= n p^2: no violation (sparse case)")
