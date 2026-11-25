from load import *
import scipy.linalg as sl
for f in FILES:
    N0,A,Bm=load(f)
    for k in (1,2,3):
        A2=[np.kron(np.eye(k,dtype=int),a) for a in A]; B2=[np.kron(np.eye(k,dtype=int),b) for b in Bm]
        N=N0*k; S,K=mats(N,A2,B2); p=np.linalg.matrix_rank((S-K).astype(float))
        print(f.split('/')[-1],'k',k,'N',N,'p',p,'N^2/2-N/2*(N0^2-2p0)/N0^2... N(N-1)/2=',N*(N-1)//2)
