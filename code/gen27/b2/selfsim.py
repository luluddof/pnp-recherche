import numpy as np, sys
src=open('dp_alpha2.py').read().split('mus=')[0]; sys.argv=[sys.argv[0]]; exec(src)
P=family('cornerC4')
def eprime(ea,eb):   # exact max over hull of min(ea u, eb v): check vertices & edges of hull (P is boundary polyline)
    U,V=P[:,0],P[:,1]
    best=np.max(np.minimum(ea*U,eb*V))
    # segments: on segment between consecutive points find crossing ea u = eb v
    U0,V0,U1,V1=U[:-1],V[:-1],U[1:],V[1:]
    d=(ea*(U1-U0)-eb*(V1-V0)); s=np.where(np.abs(d)>1e-300,(eb*V0-ea*U0)/np.where(d==0,1,d),-1)
    ok=(s>=0)&(s<=1); val=ea*(U0+s*(U1-U0))
    if ok.any(): best=max(best,np.max(val[ok]))
    return best
for alpha in (0.6,0.7,0.8,0.9):
    x=1.0
    res=[]
    for lk in np.linspace(0.5,40,400):   # x_b/x_a = 10^-lk ... use w
        k=10**(-lk); xa=1/(1+k); xb=k/(1+k)
        r=eprime(xa**-alpha, xb**-alpha)   # e(x)=x^-alpha, x=1
        res.append((r,lk))
    r,lk=min(res); print(f"alpha={alpha}: best ratio e'/x^-alpha = {r:.8f} at x_b/x_a=1e-{lk:.1f}")
