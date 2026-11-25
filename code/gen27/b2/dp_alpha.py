# DP over unbalanced grid levels (Lemma U, design level, lower-order terms dropped except 2^-x_b entry term):
# e(x) = min_w max_{(u,v) in H} min(e(w x) u, e((1-w) x) v) [+ 2^{-min(xa,xb)}]
import numpy as np, sys
def boundary(name):
    t=np.linspace(0,1,4001)
    if name=='circle':
        ang=t*np.pi/2; return np.c_[np.cos(ang),np.sin(ang)]
    if name=='line': return np.c_[t,1-t]
    if name.startswith('corner'):  # pieces (0,1),(1-1/c,1),(1,0): T<=c(2-theta)
        c=float(name[6:]); a=1-1/c
        P=np.r_[np.c_[t*a,np.ones_like(t)], np.c_[a+t*(1-a), 1-t]]
        return P
    if name.startswith('ellp'):  # l_p sphere
        p=float(name[4:]); ang=t*np.pi/2; c=np.cos(ang); s=np.sin(ang)
        r=(c**p+s**p)**(-1/p); return np.c_[r*c,r*s]
def run(name, x0=4.0, r=1.01, kmax=3000, extra=True):
    B=boundary(name)
    xs=x0*r**np.arange(kmax); le=np.zeros(kmax); le[:]=np.log(0.5)
    lx=np.log(xs)
    ws=np.r_[np.linspace(0.5,0.999,200), 1-np.logspace(-3,-8,40)]
    for k in range(1,kmax):
        x=xs[k]; best=np.log(0.5)
        for w in ws:
            xa=w*x; xb=(1-w)*x
            if xb<x0: continue
            ea=np.exp(np.interp(np.log(xa),lx[:k],le[:k])); eb=np.exp(np.interp(np.log(xb),lx[:k],le[:k]))
            # both orientations (symmetric pieces assumed? use B and swapped)
            for (U,V) in ((B[:,0],B[:,1]),(B[:,1],B[:,0])):
                val=np.max(np.minimum(ea*U,eb*V))
                if extra: val+=2.0**(-xb)
                if np.log(val)<best: best=np.log(val)
        le[k]=best
    # fit slope of -log e vs log x on last half
    h=kmax//2; sl=np.polyfit(lx[h:],-le[h:],1)[0]
    return xs[-1], np.exp(le[-1]), sl
for name in sys.argv[1:]:
    X,E,sl=run(name)
    print(f"{name:10s} x_max={X:.3g} e={E:.4g} fitted alpha={sl:.3f}", flush=True)
