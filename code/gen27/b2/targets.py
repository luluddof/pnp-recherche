# For target alpha, and weighting mu, the frontier must satisfy mu*theta + (1-mu)*T <= 2||(mu,1-mu)||_q, q=1/(1-alpha).
# Report where App B's circle theta^2+T^2=4 violates it, and how much T must drop there.
import numpy as np
th=np.linspace(0,2,200001); Tc=np.sqrt(4-th**2)
for alpha in (0.55,0.585,0.65,0.75,0.9):
    q=1/(1-alpha); print(f"alpha={alpha} (rho-equivalent {2**-alpha:.4f}; balanced needs sup(theta+T) <= {4*2**-alpha:.4f})")
    for mu in (0.5,0.6,0.7,0.8,0.9,0.95,0.99):
        nq=(mu**q+(1-mu)**q)**(1/q)
        Tmax=(2*nq-mu*th)/(1-mu)
        bad=Tc>Tmax+1e-12
        if not bad.any(): print(f"   mu={mu}: circle OK"); continue
        lo,hi=th[bad].min(),th[bad].max()
        ratio=np.max(np.where(bad, Tc/np.maximum(Tmax,1e-9),0))
        # value at the corner
        T2=(2*nq-2*mu)/(1-mu)
        print(f"   mu={mu}: violate on theta in [{lo:.3f},{hi:.3f}]  max T_circle/T_allowed={ratio:.3f}  allowed T at theta=2: {T2:.4f}")
