# Design-level DP for unbalanced grid levels (G25-P9 Lemma U), exact hull via support function:
# e' = max_{(u,v) in conv P (down-closed)} min(ea u, eb v) = min_mu h(mu) / (mu/ea + (1-mu)/eb)
import numpy as np, sys
def family(name):
    t=np.linspace(0,1,2001)
    arc=lambda a0,a1: np.c_[np.cos(a0+(a1-a0)*t), np.sin(a0+(a1-a0)*t)]
    if name=='circle': return arc(0,np.pi/2)
    if name=='line': return np.c_[t,1-t]
    if name.startswith('cornerC'):   # circle, but linear corner v=c(1-u) near (1,0) (and symmetric near (0,1))
        c=float(name[7:]); v0=2*c/(c*c+1); u0=1-v0/c
        a0=np.arctan2(v0,u0)
        seg=np.c_[1-t*(1-u0), t*v0]
        P=np.r_[seg, arc(a0,np.pi/2-a0), seg[:,::-1]]
        return P
    if name.startswith('powC'):      # circle with power corner v = k (1-u)^g near (1,0) (one corner only)
        g=float(name[4:]); 
        # junction: find s.t. curve v=k t^g meets circle; choose k so it meets at t=0.1
        tj=0.1; uj=1-tj; vj=np.sqrt(1-uj**2); k=vj/tj**g
        tt=t*tj; corner=np.c_[1-tt, k*tt**g]
        a0=np.arctan2(vj,uj); return np.r_[corner, arc(a0,np.pi/2)]
    if name.startswith('ellp'):
        p=float(name[4:]); ang=t*np.pi/2; c=np.cos(ang); s=np.sin(ang)
        r=(c**p+s**p)**(-1/p); return np.c_[r*c,r*s]
mus=np.r_[np.linspace(0,1,2001), 1-np.logspace(-4,-12,60), np.logspace(-4,-12,60)]
def run(name, x0=4.0, r=1.02, kmax=2500):
    P=family(name)
    H=np.max(np.outer(mus,P[:,0])+np.outer(1-mus,P[:,1]),axis=1)
    xs=x0*r**np.arange(kmax); le=np.full(kmax,np.log(0.5)); lx=np.log(xs)
    ws=np.r_[np.linspace(0.5,0.999,120), 1-np.logspace(-3,-10,50)]
    for k in range(1,kmax):
        x=xs[k]; best=np.log(0.5)
        for w in ws:
            for (xa,xb) in ((w*x,(1-w)*x),((1-w)*x,w*x)):
                if min(xa,xb)<x0: continue
                ea=np.exp(np.interp(np.log(xa),lx[:k],le[:k])); eb=np.exp(np.interp(np.log(xb),lx[:k],le[:k]))
                val=np.min(H/(mus/ea+(1-mus)/eb))+2.0**(-min(xa,xb))
                best=min(best,np.log(val))
        le[k]=best
    h=kmax//2; sl=np.polyfit(lx[h:],-le[h:],1)[0]
    rho=np.max(P.sum(1))/2
    return xs[-1],np.exp(le[-1]),sl,rho
for name in sys.argv[1:]:
    X,E,sl,rho=run(name)
    print(f"{name:12s} rho(balanced)={rho:.4f} log2(1/rho)={np.log2(1/rho):.3f}  x_max={X:.2g} e={E:.3g} fitted alpha={sl:.3f}",flush=True)
