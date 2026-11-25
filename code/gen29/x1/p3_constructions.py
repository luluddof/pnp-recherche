# Exact checks of P_3 constructions: C = I (ring bonds 1,1,1) and C = register transposition (13) (bonds 1,1,N^2).
import numpy as np, flint, itertools
def ranks(M):
    Z=flint.fmpz_mat([[int(v) for v in r] for r in M]); out=[Z.rank()]
    for pr in (1000003,998244353): out.append(flint.nmod_mat([[int(v)%pr for v in r] for r in M],pr).rank())
    return out
for N in range(2,6):
    D=N**3; idx=lambda a:(a[0]*N+a[1])*N+a[2]
    P=np.zeros((D,D),int); I=np.eye(D,dtype=int); T=np.zeros((D,D),int)
    for a in itertools.product(range(N),repeat=3):
        P[idx(a),idx((a[2],a[0],a[1]))]=1   # b1=a3, b2=a1, b3=a2
        T[idx(a),idx((a[2],a[1],a[0]))]=1   # b1=a3, b2=a2, b3=a1
    print(f"N={N}: rank(P3-I)={ranks(P-I)} (2/3)(N^3-N)={(2*(N**3-N))//3};  rank(P3-T13)={ranks(P-T)} (N^3-N^2)/2={(N**3-N**2)//2};  bound N^3/4={N**3/4}")
