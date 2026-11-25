# Symmetric frontier on 2*theta+T=4 (theta in [4/3,2]) and theta+2T=4: the shape of the monomial-model bound p^2 q ~ N^4.
import numpy as np
t=np.linspace(0,1,200001)
th1=4/3+t*(2/3); T1=4-2*th1           # theta in [4/3,2]
P=np.r_[np.c_[th1,T1], np.c_[T1,th1]]/2.0   # normalized (u,v)
print("balanced rho =", P.sum(1).max()/2, " -> alpha_bal =", np.log2(1/(P.sum(1).max()/2)))
for alpha in (0.6,0.8,0.9,0.95,0.99):
    q=1/(1-alpha); ok=[]
    for mu in np.linspace(0.5,0.99,491):
        h=np.max(mu*P[:,0]+(1-mu)*P[:,1]); nq=(mu**q+(1-mu)**q)**(1/q)
        if h<=nq+1e-12: ok.append(round(mu,3))
    print(f"alpha={alpha}: weightings mu that work: {ok[:3]}...{ok[-3:] if ok else ''} ({len(ok)} of 491)")
# self-similar ratio check at alpha=0.9 with exact segment geometry
def eprime(ea,eb,P=P):
    U,V=P[:,0],P[:,1]; return np.max(np.minimum(ea*U,eb*V))
for alpha in (0.8,0.9):
    best=min((eprime((1/(1+k))**-alpha,(k/(1+k))**-alpha),k) for k in 10.0**-np.linspace(0.5,20,200))
    print(f"alpha={alpha}: best self-similar ratio {best[0]:.9f} at x_b/x_a={best[1]:.2e}")
