import numpy as np
def H(p):
    p=np.clip(p,1e-15,1-1e-15); return -p*np.log2(p)-(1-p)*np.log2(1-p)
def F2(tau):   # GF(2) row rate of A_i, i=tau t, d=(1-tau)t
    d=1-tau; nu=np.linspace(max(d,1e-9),1,200001)
    return np.max(H(nu)+(1-nu)+nu*H((1-d/nu)/2))
def FQ(tau):
    z=np.linspace(0,min(tau,0.75),20001); return np.max(H(z)+z*np.log2(3))
for name,F in (("Q",FQ),("GF2",F2)):
    best=max(((F(t)+F(1-t))/4,t) for t in np.linspace(0.01,0.99,99))
    print(name,"asymptotic rho of one-bit tensor-power scheme = %.4f at tau=%.2f"%best)
