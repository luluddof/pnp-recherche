# [N] Does the exact (5,3)@4 point (G19b) contain a term G (x) G^{-1} with coefficient 1 (needed to embed it into all N with osr 3)?
import numpy as np
from scipy.optimize import least_squares
def E(i,j):
    M=np.zeros((4,4)); M[i-1,j-1]=1; return M
G=[-E(3,3)+E(4,2), -E(1,4)-E(2,3)+E(4,1), E(2,2)-E(3,1)]
Gp=[E(2,4)-E(3,3), E(1,4)-E(3,2)-E(4,1), -E(1,3)+E(2,2)]
# check (5,3): C = sum G_l(x,u) G'_l(y,v); X = swap; rank(X - C) in (x,y)|(u,v)
N=4; X=np.zeros((4,)*4)
for x in range(4):
    for y in range(4): X[x,y,y,x]=1
C=sum(np.einsum('xu,yv->xyuv',g,gp) for g,gp in zip(G,Gp))
print("rank X-C:", np.linalg.matrix_rank((X-C).reshape(16,16)))
R=np.array([g.ravel() for g in G]).T; Rp=np.array([g.ravel() for g in Gp]).T   # realigned: sum vec(G_l) vec(G'_l)^T
Q,_=np.linalg.qr(Rp)
def res(z):
    a=z[:3]+1j*z[3:]; Gm=sum(ai*g for ai,g in zip(a,G))
    try: Gi=np.linalg.inv(Gm)
    except: return np.ones(32)*1e3
    v=Gi.ravel(); r=v-Q@(Q.conj().T@v)
    return np.concatenate([r.real,r.imag])/max(1,np.linalg.norm(v))
best=9
rng=np.random.default_rng(0)
for t in range(300):
    s=least_squares(res,rng.standard_normal(6),xtol=1e-14,ftol=1e-14)
    if s.cost<best:
        best=s.cost; zb=s.x
print("min residual of G^{-1} off span(G'):", best)
a=zb[:3]+1j*zb[3:]; Gm=sum(ai*g for ai,g in zip(a,G)); Gi=np.linalg.inv(Gm)
b=np.linalg.lstsq(Rp,Gi.ravel(),rcond=None)[0]
# realigned operator Rr = R Rp^T (sum over l); coefficient lambda with rank(Rr - lam vec(G) vec(Gi)^T) = 2
Rr=R@Rp.T
for lam in [1.0]:
    print("rank after subtracting 1*G(x)G^-1 (realigned):", np.linalg.matrix_rank(Rr-lam*np.outer(Gm.ravel(),Gi.ravel()),tol=1e-8))
pinv=np.linalg.pinv(Rr)
print("needed lambda = 1/(b^T R^+ a):", 1/(Gi.ravel()@pinv@Gm.ravel()))
