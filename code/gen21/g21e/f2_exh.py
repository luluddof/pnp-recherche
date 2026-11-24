# Over F_2 itself (points, not closure): search A (3x4), B (4x3) rank 3 with max_{G,H in M_2(F_2)} rank A(GxH)B <= 2
import numpy as np, itertools
from ff import rank_mod, flat
p=2
mats=[np.array(m).reshape(2,2) for m in itertools.product(range(2),repeat=4)]
krons=[np.kron(G,H)%2 for G in mats for H in mats]
vecs=[np.array(v) for v in itertools.product(range(2),repeat=4) if any(v)]
found=[]
# A determined by kernel line u (rank 3), B by range hyperplane W (rank 3): A = basis of u-annihilator... use: rank A K B = 3 - dim(K W cap ker A) for K invertible
for u in vecs:
    A=np.array([w for w in vecs if (w@u)%2==0])  # rows spanning u^perp (dual), ker A = span u
    for h in vecs:
        B=np.array([w for w in vecs if (w@h)%2==0]).T  # columns span hyperplane h^perp
        s=max(rank_mod(A@K%2@B,2) for K in krons)
        if s<3: found.append((u,h,s))
print('pairs (u,h) with max rank <=2 over F_2 points:',len(found)); 
if found:
    u,h,s=found[0]; print('example u',u,'h',h,'s',s)
    A=np.array([w for w in vecs if (w@u)%2==0]); B=np.array([w for w in vecs if (w@h)%2==0]).T
    # Phi(X) = A X B (Fl1 block). c1 over F_2 points: rank-one X
    c1=max(rank_mod(A@np.outer(x,y)%2@B,2) for x in vecs for y in vecs)
    c2=max(rank_mod(A@K%2@B,2) for K in krons)
    F=np.zeros((4,4),dtype=int)
    for x in range(2):
        for y in range(2): F[x*2+y, y*2+x]=1
    RF=rank_mod(A@F@B,2); RI=rank_mod(A@B,2)
    print('c1',c1,'c2',c2,'R(F)',RF,'R(I)',RI,'N^2 c1 c2',4*c1*c2)
    # which pairs (p,q) this excludes for F over F_2 vs truth
