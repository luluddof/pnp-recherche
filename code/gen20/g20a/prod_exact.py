# exact check of the D=3 product point:  X_3 = YY(03,12) + [I_03 (x) F_12] + [F_03 (x) Y_12]
#   A_0 = Y_03 Y_12  (cut 0 rank (N(N-1)/2)^2),  A_1 = I_03 F_12 (cut 1 rank N^2),  A_2 = F_03 Y_12 (cut 2 rank N^2)
import numpy as np, flint, sys
from trilib import box3
for N in [2,3,4]:
    b=box3([N]*4); I=np.eye(N,dtype=np.int64)
    F=np.einsum('ad,bc->abcd',I,I); Id=np.einsum('ac,bd->abcd',I,I); Y=Id-F
    L='abcdefgh'; s1=L[0]+L[3]+L[4]+L[7]; s2=L[1]+L[2]+L[5]+L[6]
    mk=lambda O1,O2: np.einsum(f'{s1},{s2}->{L}',O1,O2)
    A=[mk(Y,Y), mk(Id,F), mk(F,Y)]
    X=b.X.astype(np.int64)
    assert (A[0]+A[1]+A[2]==X).all()   # exact integer identity, entry by entry
    rk=[flint.fmpz_mat(b.flat(A[c],c).astype(int).tolist()).rank() for c in range(3)]
    full=[[flint.fmpz_mat(b.flat(A[i],c).astype(int).tolist()).rank() for c in range(3)] for i in range(3)]
    print('N',N,'exact triple',rk,'predicted',((N*(N-1)//2)**2,N*N,N*N),'all cut ranks of the 3 terms',full,flush=True)
