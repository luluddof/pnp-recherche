# Family d: tournament / half-of-pairs matching shape.  G_l signed partial permutations with disjoint supports covering
# S = {(a,c): a<c} (or Paley tournament), G'_l = eps G_l^T.  Variants:
#  xor  : N = 2^n, M_x = {(a, a^x) : a < a^x}, signs from CD table / 1 / Clifford, x = 1..N-1 (q = N-1)
#  shift: M_l = {(a, a+l) : a+l < N}, l = 1..N-1 (q = N-1)
#  paley: r prime = 3 mod 4, M_l = {(a, a+l mod r)} for l in QR (q = (r-1)/2), full shifts
#  +diag: optionally add the diagonal term G_0 = I (G'_0 = c I)
import sys, csv, time; sys.path.insert(0, '.')
import numpy as np, scipy.sparse as sp
from ev import *
from fams import matching_from_partition, chi
from cd import table
from scan import evaluate
from sympy import isprime
def inst(N):
    if N & (N - 1) == 0 and N >= 2:
        n = N.bit_length() - 1; S, K = table(n)
        el = np.array([[(v >> k) & 1 for k in range(n)] for v in range(N)])
        B = np.tril(np.ones((n, n), int), -1); Cl = (-1) ** ((el @ B @ el.T) % 2)
        for lab, sg in (('xor-1', np.ones((N, N), int)), ('xor-CD', S), ('xor-cliff', Cl)):
            parts = [[(a, a ^ x) for a in range(N) if a < (a ^ x)] for x in range(1, N)]
            signs = [[sg[a, c] for a, c in M] for M in parts]
            for eps in (1, -1):
                yield f'{lab}-e{eps}', matching_from_partition(N, parts, signs, [eps] * len(parts))
    parts = [[(a, a + l) for a in range(N) if a + l < N] for l in range(1, N)]
    for eps in (1, -1): yield f'shift-e{eps}', matching_from_partition(N, parts, None, [eps] * len(parts))
    # both halves (all off-diagonal pairs) by shifts, eps alternating
    if isprime(N) and N % 4 == 3:
        c = chi(N); parts = [[(a, (a + l) % N) for a in range(N)] for l in range(1, N) if c[l] == 1]
        for eps in (1, -1): yield f'paley-e{eps}', matching_from_partition(N, parts, None, [eps] * len(parts))
if __name__ == '__main__':
    w = csv.writer(open('scan_famD.csv', 'a'))
    for N in map(int, sys.argv[1].split(',')):
        for lab, terms in inst(N):
            t = time.time(); p, sc = evaluate(N, terms)
            print('famD', lab, N, 'q', len(terms), 'p', p, 'scale', sc, 'N^2', N * N, '%.1fs' % (time.time() - t), flush=True)
            w.writerow(['famD', lab, N, len(terms), p, sc])
