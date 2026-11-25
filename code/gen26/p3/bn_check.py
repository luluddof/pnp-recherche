# BN model check: G_l, G'_l orthogonal (norm 1). Odd traces of Theta_s must satisfy |tr Theta_s^j| <= q^j N.
import numpy as np
from scipy.stats import ortho_group
rng=np.random.default_rng(0)
for N,q in [(5,2),(6,3),(7,4)]:
    S=np.zeros((N*N,N*N))
    for x in range(N):
        for y in range(N): S[x*N+y,y*N+x]=1
    worst=0
    for t in range(20):
        G=[ortho_group.rvs(N,random_state=rng) for _ in range(q)]; H=[ortho_group.rvs(N,random_state=rng) for _ in range(q)]
        # make some instances degenerate to get big traces: G'_l = G_l^{-1}
        if t%2: H=[g.T for g in G]
        K=sum(np.kron(g,h) for g,h in zip(G,H)); Th=K@S; Ts=(Th+Th.T)/2
        for j in (1,3,5,7):
            r=abs(np.trace(np.linalg.matrix_power(Ts,j)))/(q**j*N); worst=max(worst,r)
        # also check Theta^odd = sum (X (x) Y) S identity: tr(Theta) = sum tr(G_l H_l)
        assert abs(np.trace(Th)-sum(np.trace(g@h) for g,h in zip(G,H)))<1e-8
    print('N',N,'q',q,'max |tr Ts^j|/(q^j N) over odd j<=7 =',round(worst,4),'(must be <= 1)')
