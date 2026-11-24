import numpy as np, sys
from d2 import *
N,p,q,R = map(int, sys.argv[1:5]); fl=int(sys.argv[6]) if len(sys.argv)>6 else 1
rng=np.random.default_rng(int(sys.argv[5]))
T=swap(N); n=N*N; S=T.reshape(n,n)
np.set_printoptions(precision=3, suppress=True, linewidth=160)
for r in range(R):
    v,x,E,A,B = solve(T,N,p,q,rng.standard_normal(4*(p+q)*N*N)*0.7,4000)
    if fl==2: E=E.transpose(0,2,1,3); FF=T.transpose(0,2,1,3).reshape(n,n)
    else: FF=S
    P = FF.conj().T@E.reshape(n,n); P=(P+P.conj().T)/2
    w,V = np.linalg.eigh(P); K=V[:, w>0.5]
    print('d2',round(v,7),'coef',round(np.abs(x).max(),1),'eig P',np.round(w[w>1e-3],4))
    # symmetry: overlap of K with sym subspace
    Psym=(np.eye(n)+S)/2
    print('  dim K',K.shape[1],' tr(P_K P_sym)=',round(np.trace(K.conj().T@Psym@K).real,4))
    # Schmidt spectra of vectors in K: generic combination
    for i in range(K.shape[1]):
        print('   k%d schmidt sv'%i, np.linalg.svd(K[:,i].reshape(N,N),compute_uv=False))
    # reduced density of P_K on first and second factor
    PK=(K@K.conj().T).reshape(N,N,N,N)
    print('  tr_2 P_K eig', np.linalg.eigvalsh(np.einsum('xyuy->xu',PK)), ' tr_1 eig', np.linalg.eigvalsh(np.einsum('xyxv->yv',PK)))
