# Block (Majorana) structure: in a 2m-dim block (m qubits), a maximal anticommuting set A (2m+1 Paulis, B pairwise 1).
# g_A(w) = sum_{s in A} sig_s (-1)^{Q(s)+B(s,w)}.  Joint law of (beta_U(w), g_A(w)).
import numpy as np, itertools
from collections import Counter
def jw(m):
    # Jordan-Wigner Majoranas: gamma_{2j}=Z..Z X_j, gamma_{2j+1}=Z..Z Y_j, plus product (= +-Z..Z) => 2m+1 anticommuting
    out=[]
    for j in range(m):
        for P in 'XY':
            a=b=0
            for i in range(j): b|=1<<i
            if P=='X': a|=1<<j
            else: a|=1<<j; b|=1<<j
            out.append(a|(b<<m))
    tot=0
    for v in out: tot^=v
    out.append(tot); return out
for m in [1,2,3,4,5]:
    A=jw(m); n=4**m; msk=(1<<m)-1
    Q=lambda v: bin((v&msk)&(v>>m)).count('1')&1
    Bf=lambda u,v:(bin((u&msk)&(v>>m)).count('1')+bin((v&msk)&(u>>m)).count('1'))&1
    assert all(Bf(u,v)==1 for u,v in itertools.combinations(A,2))
    joint=Counter()
    for w in range(n):
        g=sum((-1)**(Q(s)+Bf(s,w)) for s in A); joint[(1-2*Q(w),g)]+=1
    print(m, len(A), sorted(joint.items()))
