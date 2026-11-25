# Theorem B kernel: g(sin phi)=G(phi)/G(pi/2), G(phi)=int_0^phi cos(psi) F_K(2psi) dpsi. Build exactly in sin-basis,
# convert to monomials via Chebyshev, check odd, |g|<=1, monotone, a=g(1/Q) >= 16(K+1)/(pi^4 Q) if K+1<=pi Q/2, L bound.
import numpy as np
from numpy.polynomial import chebyshev as C, polynomial as P
from fractions import Fraction
import math
def kernel(K):
    # cos psi * F_K(2psi), F_K(x)=sum_{|k|<=K}(1-|k|/(K+1)) e^{ikx}; integrate termwise -> sum_t c_t sin(t phi), t odd
    c={}
    for k in range(0,K+1):
        w=1.0 if k==0 else 2*(1-k/(K+1))
        # cos psi cos 2k psi = (cos(2k+1)psi + cos(2k-1)psi)/2 ; integral of cos(t psi) = sin(t phi)/t
        for t in (2*k+1,abs(2*k-1)):
            sgn=1.0 if (2*k-1)>=0 or t==2*k+1 else 1.0   # cos is even, so cos((2k-1)psi)=cos(|2k-1|psi)
            c[t]=c.get(t,0)+w/2/t
    # sin(t phi) with y=sin phi, t odd: sin(t phi) = (-1)^((t-1)/2) T_t(y)
    cheb=np.zeros(2*K+2)
    for t,v in c.items(): cheb[t]+=v*(-1)**((t-1)//2)
    G_end=sum(v*math.sin(t*math.pi/2) for t,v in c.items())
    cheb/=G_end
    return cheb
for K in [2,5,10,15,20]:
    ch=kernel(K); mono=C.cheb2poly(ch); L=np.sum(np.abs(mono))
    ys=np.linspace(-1,1,20001); g=C.chebval(ys,ch)
    odd=np.allclose(mono[0::2],0,atol=1e-9); mon=np.all(np.diff(g)>=-1e-12)
    Lb=2*math.pi*(1+math.log(2*K+2))*(1+math.sqrt(2))**(2*K+1)
    out=[]
    for Q in [int(math.ceil(2*(K+1)/math.pi)),20,50,100]:
        a=C.chebval(1/Q,ch); okv=(K+1)<=math.pi*Q/2
        out.append((Q,round(a*Q/(K+1),3),okv and a>=16*(K+1)/(math.pi**4*Q)))
    print(f'K={K} deg={2*K+1} odd={odd} |g|max={np.abs(g).max():.6f} g(1)={g[-1]:.6f} monotone={mon} L={L:.3e} bound={Lb:.3e} L<=bound={L<=Lb}  (Q, aQ/(K+1), a>=claim): {out}')
# asymptotic constant: K+1=0.28 ln N
for N in [1e4,1e8,1e16]:
    K=int(0.28*math.log(N))-1; ch=kernel(K); L=np.sum(np.abs(C.cheb2poly(ch)))
    print(f'N={N:.0e} K={K} L={L:.3e} sqrtN={math.sqrt(N):.3e} L*N/N^2*ln... L/sqrtN={L/math.sqrt(N):.3f}')
