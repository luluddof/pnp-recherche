# exact Bell-model optimum at N=4 (k=2): min #{phi != beta} over q-Fourier-sparse phi on F_2^4
import numpy as np, itertools
k=2; n=16
pts=np.array(list(itertools.product([0,1],repeat=4)))   # (c1,c2,d1,d2)
beta=(-1.0)**((pts[:,0]*pts[:,2]+pts[:,1]*pts[:,3])%2)
CH=(-1.0)**((pts@pts.T)%2)   # CH[xi,x] = chi_xi(x)
def B(x,y): return (x[0]*y[2]+x[1]*y[3]+y[0]*x[2]+y[1]*x[3])%2
basis=np.eye(4,dtype=int)
Sp=[]
for bits in range(1<<16):
    A=np.array([(bits>>i)&1 for i in range(16)]).reshape(4,4)
    if all(B(A@basis[i]%2,A@basis[j]%2)==B(basis[i],basis[j]) for i in range(4) for j in range(4)):
        if round(abs(np.linalg.det(A)))%2==1: Sp.append(A)
print("Sp size",len(Sp))
idx={tuple(p):i for i,p in enumerate(pts)}
perms=[]
for A in Sp:
    for l in pts:
        perms.append([idx[tuple((A.T@p+l)%2)] for p in pts])
perms=np.array(perms); print("group",len(perms))
masks=np.arange(1<<16,dtype=np.int64)
canon=masks.copy()
bitsof=((masks[:,None]>>np.arange(16))&1)
for g in perms:
    img=(bitsof<<g[None,:]).sum(1)
    canon=np.minimum(canon,img)
reps=np.unique(canon); print("orbits",len(reps))
best={}
for r in reps:
    S=[i for i in range(16) if (r>>i)&1]; q=len(S)
    if q==0: best[0]=(16,()); continue
    M=CH[S,:]  # q x 16
    bd=16; bk=None
    for Z in itertools.combinations(range(16),q):
        A=M[:,Z].T  # q x q : A[z,s]=chi_s(z)
        if abs(np.linalg.det(A))<1e-9: continue
        kap=np.linalg.solve(A,beta[list(Z)])
        phi=kap@M
        d=int((np.abs(phi-beta)>1e-9).sum())
        if d<bd: bd=d; bk=(S,kap)
    if q not in best or bd<best[q][0]: best[q]=(bd,bk)
for q in sorted(best): print("q=%2d  min p=%d"%(q,best[q][0]), best[q][1][0] if best[q][1] else "")
import pickle; pickle.dump(best,open("bell_n4.pkl","wb"))
