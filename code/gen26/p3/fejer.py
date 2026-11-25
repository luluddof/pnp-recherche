# Theorem B kernel: g(sin phi)=G(phi)/G(pi/2), G(phi)=int_0^phi cos(psi) F_K(2 psi) dpsi.  Check a=g(1/Q) and L.
import numpy as np
from numpy.polynomial import chebyshev as C
def build(K):
    # cos(psi)*F_K(2psi) = sum_{|j|<=K} w_j cos(2 j psi) cos(psi), w_j = 1-|j|/(K+1)
    # = sum_j w_j [cos((2j+1)psi)+cos((2j-1)psi)]/2 ; collect coefficient a_t of cos(t psi), t odd>0
    a={}
    for j in range(-K,K+1):
        w=1-abs(j)/(K+1)
        for t in (2*j+1,2*j-1):
            a[abs(t)]=a.get(abs(t),0)+w/2
    # G(phi)=sum a_t sin(t phi)/t ; sin(t phi)=(-1)^((t-1)/2) T_t(sin phi)
    deg=2*K+1; beta=np.zeros(deg+1)
    for t,v in a.items(): beta[t]+=v/t*(-1)**((t-1)//2)
    Gmax=C.chebval(1.0,beta); beta/=Gmax
    return beta
for K,Q in [(5,20),(10,50),(20,100),(20,40)]:
    beta=build(K); ys=np.linspace(-1,1,20001); gv=C.chebval(ys,beta)
    mono=C.cheb2poly(beta); L=np.sum(np.abs(mono))
    a=C.chebval(1/Q,beta); mon=np.all(np.diff(gv)>=-1e-12)
    print('K',K,'Q',Q,'max|g|=%.4f'%np.max(np.abs(gv)),'monotone',mon,'a=g(1/Q)=%.4f'%a,'claimed >= %.4f'%(16*(K+1)/(np.pi**4*Q)) if K+1<=np.pi*Q/2 else '','L=%.3g'%L,'claimed <= %.3g'%(2*np.pi*(1+np.log(2*K+2))*(1+np.sqrt(2))**(2*K+1)))
