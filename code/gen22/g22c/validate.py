import json, sys, time, numpy as np
sys.path.insert(0, '.')
from ev import *
from cd import table
B = '/home/user/pnp-recherche/code/'
for f in ['gen19/g19b/exact_4_5_3.json', 'gen20/g20b/exact_m_4_4_4.json', 'gen19/g19b/exact_5_6_6.json',
          'gen19/g19b/exact_5_7_4.json', 'gen19/g19b/exact_6_7_7.json']:
    J = json.load(open(B + f)); N = J['N']
    terms = [(np.array(g, dtype=np.int64), np.array(gp, dtype=np.int64)) for g, gp in zip(J['G'], J['Gp'])]
    print(f, 'N', N, 'claimed', (J['p'], J['q']), 'eval', p_of(N, terms), 'dense', dense_check(N, terms), 'q', len(terms))

def cd_terms(n):
    S, K = table(n); N = 1 << n
    T = []
    for x in range(N):
        G = np.zeros((N, N), dtype=np.int64); Gp = np.zeros((N, N), dtype=np.int64)
        for a in range(N):
            c = a ^ x; G[a, c] = S[a, c]          # G_x(a,c) = s(a,c)[x = a^c]
        for b in range(N):
            dd = b ^ x; Gp[b, dd] = S[dd, b]      # G'_x(b,d) = s(d,b)[x = d^b]  (times 1/2 -> d=2)
        T.append((G, Gp))
    return T
for n in range(5):
    N = 1 << n; T = cd_terms(n)
    t = time.time(); r = p_of(N, T, d=2)
    print('CD N', N, 'q', N, 'p', r, 'dense', dense_check(N, T, d=2) if N <= 16 else '-', '%.2fs' % (time.time() - t))
