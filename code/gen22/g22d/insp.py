import numpy as np, sys
from d2 import *
N,p,q,R = map(int, sys.argv[1:5]); rng=np.random.default_rng(int(sys.argv[5]) if len(sys.argv)>5 else 0)
T=swap(N)
for r in range(R):
    v,x,E,A,B = solve(T,N,p,q,rng.standard_normal(4*(p+q)*N*N)*0.7,4000)
    a,b,G,H = unpack(x,N,p,q)
    f=lambda M: np.round(np.linalg.svd(M,compute_uv=False),4)
    print('d2',round(v,8),'|A|^2',round(np.linalg.norm(A)**2,4),'|B|^2',round(np.linalg.norm(B)**2,4),'<A,B>',np.round(np.vdot(A,B),4),'<F,A>',np.round(np.vdot(A,T),4),'<F,B>',np.round(np.vdot(B,T),4))
    for k in range(p): print('  a sv',f(a[k]),' b sv',f(b[k]))
    for k in range(q): print('  G sv',f(G[k]),' H sv',f(H[k]))
    print('  E Fl1',f(E.reshape(N*N,N*N))[:6],'Fl2',f(E.transpose(0,2,1,3).reshape(N*N,N*N))[:6],'Fl3',f(E.transpose(0,3,1,2).reshape(N*N,N*N))[:6])
    print('  A+B Fl1',f((A+B).reshape(N*N,N*N)), ' Fl3', f((A+B).transpose(0,3,1,2).reshape(N*N,N*N)))
