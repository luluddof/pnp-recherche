# cell-graph block structure of Theta = K S on the N^2 cells: connected components of supp(Theta) (symmetrised),
# per block: size c, rank of (1-Theta) on the block, eigenvalues, whether block = I - D J D (rank-1 defect).
import sys, numpy as np
sys.path.insert(0,'/home/user/pnp-recherche/code/gen26/p3')
from load import load, mats, FILES
from collections import Counter
from scipy.sparse.csgraph import connected_components
fs=sys.argv[1:] or FILES+['exact_6_8_6.json']
for f in fs:
    N,A,B=load(f); S,K=mats(N,A,B); T=(K@S).astype(float); n=N*N
    adj=((T!=0)|(T.T!=0)).astype(int); nc,lab=connected_components(adj,directed=False)
    stats=Counter()
    for c in range(nc):
        idx=np.where(lab==c)[0]; Bk=T[np.ix_(idx,idx)]; r=np.linalg.matrix_rank(np.eye(len(idx))-Bk)
        ev=np.round(np.sort(np.linalg.eigvals(Bk).real),3)
        stats[(len(idx),r,tuple(ev))]+=1
    tot=sum(k[1]*v for k,v in stats.items())
    print(f.split('/')[-1],"N",N,"q",len(A),"p",tot,"blocks (size, rank(1-Theta_block), eig): count")
    for k,v in sorted(stats.items()): print("   ",k[0],k[1],k[2] if k[0]<=6 else str(k[2][:4])+'...',"x",v)
