# remove diagonal + per-term scaling gauge by log-linear least squares so that nonzero magnitudes become 1, then round and save.
import numpy as np, sys
f,N,q,out=sys.argv[1],int(sys.argv[2]),int(sys.argv[3]),sys.argv[4]; thr=float(sys.argv[5]) if len(sys.argv)>5 else 0.02
z=np.load(f); s=N*N; G=z[:q*s].reshape(q,N,N).copy(); H=z[q*s:].reshape(q,N,N).copy()
G[abs(G)<thr]=0; H[abs(H)<thr]=0
rows=[];rhs=[]; nv=2*N+q   # xi(N), eta(N), tau(q)
for l in range(q):
    for a,b in zip(*np.nonzero(G[l])):
        r=np.zeros(nv); r[2*N+l]=1; r[a]+=1; r[N+b]-=1; rows.append(r); rhs.append(-np.log(abs(G[l,a,b])))
    for a,b in zip(*np.nonzero(H[l])):
        r=np.zeros(nv); r[2*N+l]=-1; r[N+a]+=1; r[b]-=1; rows.append(r); rhs.append(-np.log(abs(H[l,a,b])))
A=np.array(rows); w,res,rk,_=np.linalg.lstsq(A,np.array(rhs),rcond=None)
xi,eta,tau=w[:N],w[N:2*N],w[2*N:]
X=np.diag(np.exp(xi)); Y=np.diag(np.exp(eta))
G2=np.array([np.exp(tau[l])*X@G[l]@np.linalg.inv(Y) for l in range(q)]); H2=np.array([np.exp(-tau[l])*Y@H[l]@np.linalg.inv(X) for l in range(q)])
print("residual of log fit:", float(np.abs(A@w-np.array(rhs)).max()))
print("magnitudes after:", sorted(set(np.round(np.abs(np.concatenate([G2[G2!=0],H2[H2!=0]])),3))))
np.save(out,np.concatenate([G2.ravel(),H2.ravel()]))
