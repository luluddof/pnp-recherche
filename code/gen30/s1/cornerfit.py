# best known Full (= Bell, Majorana / exact) points; normalized N-units: q = N^eta, p = N^{2-delta}
import numpy as np
pts = {2: [(1,6),(4,4),(6,1)], 3: [(1,28),(4,24),(6,19),(8,8),(28,1)], 4: [(1,120),(10,46),(46,10)],
       5: [(12,232),(67,67),(232,12)], 6: [(14,1093),(92,378),(378,92)], 7: [(16,4944),(121,1941),(576,576),(1941,121)]}
print(' k   q      p    eta   delta  delta/eta  delta/eta^2   line(delta=eta) lower bd on p: N^2/(4q)')
corner = []
for k, L in pts.items():
    N = 2**k
    for q, p in L:
        eta = np.log(q)/np.log(N); dl = 2 - np.log(p)/np.log(N)
        if q > 1: print(f'{k:2d} {q:4d} {p:6d}  {eta:.3f}  {dl:.3f}   {dl/eta:.3f}     {dl/eta**2:.3f}        {N*N/(4*q):.1f}')
        if q == 2*k + 2 or (k <= 3 and q == 2*k+2): corner.append((k, eta, dl))
c = np.array([(e, d) for _, e, d in corner])
print('corner points (q = 2k+2, Majorana d=1):', [(k, round(e,3), round(d,3)) for k, e, d in corner])
# fits delta = a eta^b
b, la = np.polyfit(np.log(c[:,0]), np.log(c[:,1]), 1); print(f'power fit on corner points: delta ~ {np.exp(la):.3f} eta^{b:.2f}')
# asymptotic Majorana corner from H1: delta ~ ln2 eta^2 / ln^2(1/eta)
for eta in (0.3, 0.1, 0.03, 0.01):
    print(f'  H1 asymptotic: eta={eta}: delta ~ {np.log(2)*eta**2/np.log(1/eta)**2:.2e}  (line: {eta})')
