import numpy as np, sys
src=open('dp_alpha2.py').read().split('for name in sys.argv')[0]
sys_argv=sys.argv; sys.argv=[sys.argv[0]]; exec(src)
def run2(name, x0=4.0, r=1.05, kmax=3000):
    P=family(name)
    H=np.max(np.outer(mus,P[:,0])+np.outer(1-mus,P[:,1]),axis=1)
    xs=x0*r**np.arange(kmax); le=np.full(kmax,np.log(0.5)); lx=np.log(xs)
    ws=np.r_[np.linspace(0.5,0.999,60), 1-np.logspace(-3,-25,90)]
    for k in range(1,kmax):
        x=xs[k]; best=np.log(0.5)
        for w in ws:
            xa=w*x; xb=(1-w)*x
            if xb<x0: continue
            ea=np.exp(np.interp(np.log(xa),lx[:k],le[:k])); eb=np.exp(np.interp(np.log(xb),lx[:k],le[:k]))
            for (A,B) in ((ea,eb),(eb,ea)):
                val=np.min(H/(mus/A+(1-mus)/B))+2.0**(-xb)
                best=min(best,np.log(val))
        le[k]=best
    for frac in (0.25,0.5,0.75,0.9):
        i=int(kmax*frac); j=min(kmax-1,i+kmax//10)
        print(f"  {name} local alpha around x=1e{np.log10(xs[i]):.0f}: {-(le[j]-le[i])/(lx[j]-lx[i]):.3f}")
for n in sys_argv[1:]: run2(n)
