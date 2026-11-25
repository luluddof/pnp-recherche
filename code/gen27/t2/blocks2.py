import sys, numpy as np
from lib import load, swapm, Kmat
from collections import Counter
from scipy.sparse.csgraph import connected_components
for f in sys.argv[1:]:
    N,E,C=load(f); S=swapm(N); T=(Kmat(E,C)@S).astype(float); n=N*N
    adj=((T!=0)|(T.T!=0)).astype(int); nc,lab=connected_components(adj,directed=False); st=Counter()
    for c in range(nc):
        idx=np.where(lab==c)[0]; Bk=T[np.ix_(idx,idx)]; r=np.linalg.matrix_rank(np.eye(len(idx))-Bk)
        ev=np.round(np.sort(np.linalg.eigvals(Bk).real),2); st[(len(idx),r,tuple(ev) if len(idx)<=4 else 'big')]+=1
    ent=sorted(set(np.concatenate([np.ravel(x) for x in E+C]).tolist()))
    mono=all(((np.abs(x)>0).sum(0)<=1).all() and ((np.abs(x)>0).sum(1)<=1).all() for x in E+C)
    print(f.split('/')[-1],'p=',sum(k[1]*v for k,v in st.items()),'entries',ent,'MC(monomial)',mono,'symTheta',bool((T==T.T).all()),'trTheta',int(T.trace()))
    for k,v in sorted(st.items(),key=lambda t:(t[0][0],t[0][1])): print('   size',k[0],'cost',k[1],k[2],'x',v)
