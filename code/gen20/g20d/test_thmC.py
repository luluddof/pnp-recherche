# Theorem C (a),(b) as fixed linear maps; check profile (c1,c2,R) = (1,k^2,kN) and (k^2,1,kN).
import numpy as np, sys
sys.path.insert(0,'.')
from phi import *
for N in [3,4,5]:
    for k in range(1,N+1):
        I = np.eye(N, dtype=np.int64).astype(object)
        pi = np.zeros((k,N),dtype=object); sig = np.zeros((N,k),dtype=object)
        for i in range(k): pi[i,i]=1; sig[i,i]=1
        A = np.kron(pi, I); B = np.kron(I, sig)        # (pi (x) 1) X (1 (x) sigma)
        a = profile([(1,A,B)], N); b = profile([(2,A,B)], N)
        print(N,k,'(a) c1,c2,R =',a,' bound p + k^2 q >= kN;  (b):',b, ' R^2/(N^2 c1 c2)=%.3f %.3f'%(a[2]**2/(N*N*a[0]*a[1]), b[2]**2/(N*N*b[0]*b[1])))
