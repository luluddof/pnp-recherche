# Theorem A self-contained kernel: G'(psi)=cos(psi) F_K(2psi)^r ; check |g|<=1, 1-g(delta) <= u^{2r-1} pi/(1.4(2r-1)), L bound.
import numpy as np
from numpy.polynomial import chebyshev as C
def build(K,r):
    # Fourier coefficients of F_K^r via polynomial convolution of Fejer coefficients (in variable s)
    w=np.array([1-abs(j)/(K+1) for j in range(-K,K+1)])
    c=np.array([1.0])
    for _ in range(r): c=np.convolve(c,w)          # coefficients of e^{i j s}, j=-Kr..Kr
    M=K*r; a={}
    for idx,v in enumerate(c):
        j=idx-M
        for t in (2*j+1,2*j-1): a[abs(t)]=a.get(abs(t),0)+v/2
    deg=2*M+1; beta=np.zeros(deg+1)
    for t,v in a.items(): beta[t]+=v/t*(-1)**((t-1)//2)
    return beta/C.chebval(1.0,beta)
for Q,K,r in [(3,25,2),(3,25,3),(5,40,3),(5,40,4)]:
    d=1/Q; beta=build(K,r); ys=np.linspace(-1,1,40001); gv=C.chebval(ys,beta)
    u=np.pi**2/(4*(K+1)*d); bound=u**(2*r-1)*np.pi/(1.4*(2*r-1))
    k=2*K*r+1; Lb=1.82*(1+np.log(k))*(np.pi**2/4)**r*(1+np.sqrt(2))**k
    L=np.sum(np.abs(C.cheb2poly(beta))) if k<160 else float('nan')
    print('Q',Q,'K',K,'r',r,'deg',k,'max|g|=%.6f'%np.abs(gv).max(),'1-g(1/Q)=%.2e'%(1-C.chebval(d,beta)),'bound %.2e'%bound,'(u=%.3f)'%u,'L=%.3g'%L,'Lbound %.3g'%Lb)
