from load import *
for f in FILES:
    N0,A,Bm=load(f)
    for m in (1,2,3):
        for mode in ('IA','AI'):
            if mode=='AI': A2=[np.kron(a,np.eye(m,dtype=int)) for a in A]; B2=[np.kron(b,np.eye(m,dtype=int)) for b in Bm]
            else: A2=[np.kron(np.eye(m,dtype=int),a) for a in A]; B2=[np.kron(np.eye(m,dtype=int),b) for b in Bm]
            N=N0*m; S,K=mats(N,A2,B2); p=np.linalg.matrix_rank((S-K).astype(float))
            if mode=='AI' or m==1: print(f.split('/')[-1],'A(x)I_%d'%m,'N',N,'q',len(A),'p',p,'N^2/2-p=',N*N/2-p,' N(N-1)/2-p=',N*(N-1)//2-p)
