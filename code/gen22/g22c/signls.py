# Family f/a: local search over sign tables s(i,j) on Z_2^n (structure-tensor ansatz, q = N slices, scale 1/2):
# minimise p = rank(2F - sum_x T_x (x) T'_x).  Start from Cayley-Dickson (or random).  Exact ranks (blocks over Q / mod p).
import sys, time, json; sys.path.insert(0,'.')
import numpy as np, scipy.sparse as sp
from ev import *
from cd import table
n = int(sys.argv[1]); start = sys.argv[2]; T_LIM = float(sys.argv[3]); seed = int(sys.argv[4])
N = 1 << n; a = np.arange(N); rng = np.random.default_rng(seed)
S, K = table(n)
if start == 'rand':
    S = rng.choice([-1, 1], size=(N, N)); S[0, :] = 1; S[:, 0] = 1
def ev(S):
    T = []
    for x in range(N):
        c = a ^ x
        T.append((sp.csr_matrix((S[a, c], (a, c)), shape=(N, N)), sp.csr_matrix((S[c, a], (a, c)), shape=(N, N))))
    return p_of(N, T, d=2, exact_limit=0)[0]   # mod-prime ranks during search
cur = ev(S); best = cur; bestS = S.copy(); t0 = time.time(); it = 0; Temp = 2.0
while time.time() - t0 < T_LIM:
    it += 1
    i, j = rng.integers(1, N, 2)
    S[i, j] *= -1
    new = ev(S)
    if new <= cur or rng.random() < np.exp(-(new - cur) / Temp):
        cur = new
        if cur < best:
            best = cur; bestS = S.copy(); print('it', it, 'best', best, '%.0fs' % (time.time() - t0), flush=True)
    else:
        S[i, j] *= -1
    Temp = max(0.05, 2.0 * (1 - (time.time() - t0) / T_LIM))
T = []
for x in range(N):
    c = a ^ x; T.append((sp.csr_matrix((bestS[a, c], (a, c)), shape=(N, N)), sp.csr_matrix((bestS[c, a], (a, c)), shape=(N, N))))
print('FINAL', n, start, 'best', best, 'exactQ', p_of(N, T, d=2), 'iters', it)
json.dump({'n': n, 'S': bestS.tolist(), 'p': best}, open(f'signls_{n}_{start}_{seed}.json', 'w'))
