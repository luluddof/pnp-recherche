import numpy as np, sys
from d2 import *
N,p,q,R = map(int, sys.argv[1:5]); rng=np.random.default_rng(int(sys.argv[5]) if len(sys.argv)>5 else 0)
T=swap(N); n=N*N
Fm=T.reshape(n,n)
def Rl(X): return X.transpose(0,2,1,3)
for r in range(R):
    v,x,E,A,B = solve(T,N,p,q,rng.standard_normal(4*(p+q)*N*N)*0.7,4000)
    E1=E.reshape(n,n); E2=Rl(E).reshape(n,n); F2=Rl(T).reshape(n,n)
    np.set_printoptions(precision=3, suppress=True, linewidth=150)
    print('d2',round(v,6),'coef',round(np.abs(x).max(),1))
    for name,M,FF in [('Fl1',E1,Fm),('Fl2',E2,F2)]:
        print(' ',name,'eig(FE)',np.round(np.linalg.eigvals(FF.conj().T@M),3), ' herm-dev', round(np.linalg.norm(FF.conj().T@M-(FF.conj().T@M).conj().T),4))
        print('    eig(E F^*)', np.round(np.linalg.eigvals(M@FF.conj().T),3))
