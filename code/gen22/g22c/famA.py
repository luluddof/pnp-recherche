# Family a: Cayley-Dickson structure-tensor slices, q = number of slices used; p = rank(2F - sum_x T_x (x) T'_x)
import sys, time, csv; sys.path.insert(0,'.')
import numpy as np, scipy.sparse as sp
from ev import *
from cd import table
def cd_terms_sp(n, xs=None):
    S, K = table(n); N = 1 << n; a = np.arange(N)
    T = []
    for x in (range(N) if xs is None else xs):
        c = a ^ x
        G = sp.csr_matrix((S[a, c], (a, c)), shape=(N, N))
        Gp = sp.csr_matrix((S[c, a], (a, c)), shape=(N, N))   # G'_x(b,d) = s(d,b), d = b^x
        T.append((G, Gp))
    return T
if __name__ == '__main__':
  out = csv.writer(open('famA_cd.csv', 'a'))
  for n in map(int, sys.argv[1:]):
    N = 1 << n
    for q in sorted(set([N] + [N // 2, 3 * N // 4, N - 1, 1, 2, 4])):
        if q > N or q < 1: continue
        t = time.time(); r, m = p_of(N, cd_terms_sp(n, range(q)), d=2)
        print(N, q, r, m, '%.1fs' % (time.time() - t), flush=True); out.writerow(['CD', N, q, r, m])
