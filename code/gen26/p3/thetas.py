from load import *
for f in FILES:
    N,A,Bm=load(f); S,K=mats(N,A,Bm); q=len(A)
    Th=(K@S).astype(float)          # Theta = K Swap ; rank(1-Theta)=rank(Swap-K)
    p=np.linalg.matrix_rank(np.eye(N*N)-Th)
    Ts=(Th+Th.T)/2; ev=np.linalg.eigvalsh(Ts)
    n=N*N; m=n-p
    nge1=int(np.sum(ev>=1-1e-9)); nle_1=int(np.sum(ev<=-1+1e-9))
    tr=[round(np.trace(np.linalg.matrix_power(Ts,j)),3) for j in (1,3,5,7)]
    print(f.split('/')[-1],'N',N,'q',q,'p',p,'m',m,'#ev>=1',nge1,'#ev<=-1',nle_1,'||Ts||=%.3f'%max(abs(ev)),'odd traces',tr,'bounds q^jN',[q**j*N for j in (1,3,5,7)])
