import sys, math, numpy as np, importlib.util
spec=importlib.util.spec_from_file_location('appb','/home/user/pnp-recherche/code/gen20/g20b/appb_exact.py'); appb=importlib.util.module_from_spec(spec); spec.loader.exec_module(appb)
from collections import Counter
h=int(sys.argv[1]); N,mus,I,T=appb.pieces(h); d=math.factorial(h); n=N*N
def mat(T4): # rows (x,y) cols (u,v)
    return np.array([[T4[x][y][u][v] for u in range(N) for v in range(N)] for x in range(N) for y in range(N)],dtype=float)
S=np.zeros((n,n))
for x in range(N):
    for y in range(N): S[x*N+y,y*N+x]=1
A={s:mat(T[s]) for s in I}
def rc(M): return np.linalg.matrix_rank(M.reshape(N,N,N,N).transpose(0,2,1,3).reshape(n,n),tol=1e-7*max(1,abs(M).max()))
R=np.zeros((n,n))
for s0 in I+[I[-1]+1]:
    C=d*S-R; K=C/d; Th=K@S; p=np.linalg.matrix_rank(np.eye(n)-Th,tol=1e-9); q=rc(C)
    if 0<p<n and 0<q<n:
        Ts=(Th+Th.T)/2; ev=np.linalg.eigvalsh(Ts)
        # Q: nuclear-type norm: min over decompositions is >= ; use SVD of realigned K: sum sigma_i ||U_i|| ||V_i||
        Kr=K.reshape(N,N,N,N).transpose(0,2,1,3).reshape(n,n); U,sg,Vt=np.linalg.svd(Kr)
        Qs=sum(sg[i]*np.linalg.norm(U[:,i].reshape(N,N),2)*np.linalg.norm(Vt[i].reshape(N,N),2) for i in range(q))
        print(f"h={h} N={N} s0={s0} p={p} q={q} pq/n={p*q/n:.2f} ||Ts||={max(abs(ev)):.2f} Q_svd={Qs:.2f} trTheta={np.trace(Th):.2f} #ev>=1:{int(np.sum(ev>=1-1e-9))} (m={n-p}) #ev<0:{int(np.sum(ev<-1e-9))} min ev={ev.min():.2f} hist(rest)={np.histogram(ev[abs(ev-1)>1e-8],bins=[-1e9,-10,-3,-1.5,-0.5,0.5,1.5,3,10,1e9])[0].tolist()}")
    if s0 in A: R=R+A[s0]
