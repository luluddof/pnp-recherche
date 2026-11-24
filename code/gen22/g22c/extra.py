import sys, time, csv; sys.path.insert(0, '.')
import numpy as np, scipy.sparse as sp, flint
from ev import *
from famA import cd_terms_sp
from scan import evaluate
job = sys.argv[1]
if job == 'kron':   # Kronecker square of the 8-square identity: C = C1 (x) C2 only, N = 64, q = 64
    T = cd_terms_sp(3); terms = [(sp.kron(G1, G2, format='csr'), sp.kron(H1, H2, format='csr')) for G1, H1 in T for G2, H2 in T]
    print('kron8x8 N=64 q', len(terms), 'p(scale 1/4)', p_of(64, terms, d=4), flush=True)
    print('best_scale', best_scale(64, terms, maxblock=4096), flush=True)
if job == 'cd256':
    t = time.time(); print('CD N=256 p', p_of(256, cd_terms_sp(8), d=2), '%.0fs' % (time.time() - t), flush=True)
if job == 'trunc':   # CD slices x < q (q = 2^k: subalgebra slices) and x in random subsets, best scale
    w = csv.writer(open('famA_trunc.csv', 'a'))
    for n in (3, 4, 5):
        N = 1 << n
        for q in range(1, N + 1):
            p, sc = evaluate(N, cd_terms_sp(n, range(q)))
            print('CDtrunc', N, q, p, sc, flush=True); w.writerow(['CDtrunc', N, q, p, sc])
if job == 'dense32':
    T = cd_terms_sp(5); N = 32
    t = time.time(); print('dense Q rank CD N=32:', dense_check(N, T, d=2), '%.0fs' % (time.time() - t), flush=True)
