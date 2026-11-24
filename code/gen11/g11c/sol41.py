# Find + independently verify a non-commutative J(t) for (l,d)=(4,1) with rank profile (2,3,3,3) (total 11 < 12).
import numpy as np, sys
from feas import Prob
l,d=4,1; tau=[1.0,1.7,2.9,4.3]; prof=[2,3,3,3]
P=Prob(l,d,tau,prof)
rng=np.random.default_rng(12345)
for it in range(200):
    o=P.lm(rng.normal(size=P.np_),iters=800,tol=1e-28)
    print(it,"cost %.2e"%o.cost,flush=True)
    if o.cost<1e-24: break
p=o.x
# polish with Gauss-Newton (min-norm steps)
for _ in range(5):
    r=P.res(p); J=P.jac(p); p=p-np.linalg.lstsq(J,r,rcond=None)[0]
r=P.res(p); J=P.jac(p); sv=np.linalg.svd(J,compute_uv=False)
print("residual 2-norm %.3e ; #res %d #params %d ; sigma_min(J) over first #res sing. values = %.3e, sigma_max %.3e"%(np.linalg.norm(r),P.nres,P.np_,sv[P.nres-1],sv[0]))
np.save("sol41_2333.npy",p)
# independent reconstruction: A1[x,y] from the blocks, J(t)=X+t*A1
n=1<<l; X=np.zeros((n,n)); 
for x in range(n): X[x,(n-1)^x]=1
e=P.evals(p); A1=np.zeros((n,n)); spread=0
for x in range(n):
    for y in range(n):
        D=[k for k in range(l) if ((x^y)>>k)&1]
        if not D: continue
        vals=[(e[P.eidx[(k,x,y)]]-X[x,y])/tau[k] for k in D]
        spread=max(spread,max(vals)-min(vals)); A1[x,y]=np.mean(vals)
print("max inconsistency of slopes across directions: %.3e"%spread)
np.save("sol41_2333_A1.npy",A1)
for k in range(l):
    Jt=X+tau[k]*A1
    r0=[x for x in range(n) if not (x>>k)&1]; r1=[x for x in range(n) if (x>>k)&1]
    for nm,(R,C) in (("B01",(r0,r1)),("B10",(r1,r0))):
        s=np.linalg.svd(Jt[np.ix_(R,C)],compute_uv=False)
        print("dir",k,nm,"sing vals:",np.array2string(s,precision=3,suppress_small=False,max_line_width=200))
# commutativity test: is G(t)=J(t)X a group-circulant?  and does A1 X commute with the X_j's?
N=-A1@X
dev=0
for a in range(1,n):
    Pm=np.zeros((n,n)); 
    for x in range(n): Pm[x,x^a]=1
    dev=max(dev,np.abs(Pm@N@Pm.T-N).max())
print("max deviation of N=-A1 X from Z_2^l-circulant: %.3e  (|N|max %.3e)"%(dev,np.abs(N).max()))
print("N symmetric? max|N-N^T| = %.3e"%np.abs(N-N.T).max())
