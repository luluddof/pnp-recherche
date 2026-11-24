import numpy as np, scipy.optimize as so, sys
N=int(sys.argv[1]); w=int(sys.argv[2]); rng=np.random.default_rng(2)
sz=[N*w*w,N*w*w,N*w*w,N*w*w]
def unpack(z):
    o=0;out=[]
    for shp in [(N,w,w),(N,w,w),(N,w*w),(N,w*w)]:
        n=int(np.prod(shp)); out.append(z[o:o+n].reshape(shp)); o+=n
    return out
def res(z):
    A,B,Y,U=unpack(z)
    T=np.einsum('uij,xik,vjl,ykl->xvuy',U.reshape(N,w,w),A,B,Y.reshape(N,w,w))
    E=np.einsum('xv,uy->xvuy',np.eye(N),np.eye(N))
    return np.concatenate([(T-E).ravel(),(U@Y.T-np.eye(N)).ravel()])
best=None
for t in range(60):
    s=so.least_squares(res,rng.normal(size=sum(sz)),xtol=1e-15,ftol=1e-15,gtol=1e-15)
    if best is None or s.cost<best.cost: best=s
    if s.cost<1e-26: break
print("N",N,"w",w,"restarts",t+1,"max residual",np.abs(res(best.x)).max())
