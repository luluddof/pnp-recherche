# Integer-x DP allowing tiny children (non-self-similar schedules). E[x] = best e at size <= x (running min).
# split: big child x_a = x - y - ov, small child y; ov = log t_a + log t_b (zeta=x^-1.5 -> t ~ x^4 -> ov = 8 log2 x) if ovh
import numpy as np, sys
sys.path.insert(0,'.'); from r1dp import Hull, fam
def run(name, X, ovh):
    H=None if name=='circle' else Hull(fam(name))
    E=np.full(X+1,0.5)
    Yall=np.unique(np.r_[np.arange(1,64),np.round(np.logspace(np.log10(64),np.log10(X),300)).astype(int)])
    rep=[10**k for k in range(2,int(np.log10(X))+1)]; out=[]
    for x in range(2,X+1):
        ov=8*np.log2(x) if ovh else 0.0
        Y=Yall[Yall<=x//2]
        xa=(x-Y-np.ceil(ov)).astype(int); ok=xa>=1
        if not ok.any(): E[x]=E[x-1]; continue
        Y=Y[ok]; xa=xa[ok]
        ea=E[xa]; eb=E[Y]
        if ovh:
            z=x**-1.5
            ea=ea*(1+2*z+2.0**(-Y.astype(float))); eb=eb*(1+2*z+2.0**(-xa.astype(float)))
        if H is None: F=1/np.sqrt(1/ea**2+1/eb**2)
        else: F=np.minimum(H.s_star(ea,eb),H.s_star(eb,ea))
        if ovh: F=F+2.0**(1-np.minimum(Y,xa).astype(float))
        i=np.argmin(F)
        E[x]=min(E[x-1],F[i])
        if x in rep or x%100000==0: print(f"  {name} x={x:.0e} e={E[x]:.4e} e*x={E[x]*x:.4g} best y={Y[i]}",flush=True)
    return out,E
if __name__=='__main__':
    X=int(float(sys.argv[1])); ovh=int(sys.argv[2])
    for name in sys.argv[3:]:
        out,E=run(name,X,ovh); print(name,'ovh',ovh,flush=True); print('  '+'\n  '.join(out),flush=True)
        np.save(f'E_{name}_{ovh}.npy',E)
