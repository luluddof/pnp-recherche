# min pq/N^2 over integer p,q>=1 with p >= max_k k(N-kq), q >= max_k k(N-kp) (Theorem 2 of G18g)
import numpy as np
worst=1e9
for N in range(2,201):
    k=np.arange(0,N+1)
    lp=lambda q: int((k*(N-k*q)).max())
    best=1e9
    for q in range(1,N*N+1):
        p=max(1,lp(q))
        if q< max(1,int((k*(N-k*p)).max())): continue
        best=min(best,p*q/N**2)
        if q>2*N: break
    worst=min(worst,best)
print('min over N<=200 of min pq/N^2 =',worst,' 3/16=',3/16,' 2/9=',2/9)
