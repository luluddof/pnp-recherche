# Independent DP over unbalanced Lemma-U trees. Rigorous rounding: value stored at grid x_k is achieved at some size <= x_k
# (children sizes rounded DOWN, parent size = sum (+ grouping overhead) <= x_k). e(x)=1/2 brute-force base everywhere.
# Level map geometric (ray/hull), circle exact (Pythagoras).  Option ovh=1: multiplicative (1+2zeta+2^-xmin),
# additive 2^(1-xmin) entry term, and x' = xa+xb+log t_a+log t_b with zeta=1/x', t = x'^2 ln(2^(2x')) / 2 (Lemma 8(ii)).
import numpy as np, sys
sys.path.insert(0,'.'); from r1dp import Hull, fam
def run(name, ovh=0, x0=2.0, r=1.02, xmax=1e40, report=(1e2,1e3,1e4,1e6,1e8,1e12,1e16,1e20,1e30,1e40)):
    if name!='circle': H=Hull(fam(name))
    n=int(np.log(xmax/x0)/np.log(r))+1
    xs=x0*r**np.arange(n); le=np.full(n,np.log(0.5))
    for k in range(1,n):
        xk=xs[k]
        # candidate pairs: i = index of larger child (both orders evaluated); j = largest with x_j <= budget
        i=np.arange(k)
        xa=xs[i]
        if ovh:
            # overhead: approximate log t <= 3 log2(xk)+2 per child (zeta=1/xk)
            lt=3*np.log2(xk)+2
            budget=xk-xa-2*lt
        else: budget=xk-xa
        j=np.searchsorted(xs,budget,side='right')-1
        ok=(j>=0)&(j<k)
        if not ok.any(): le[k]=min(le[k-1],np.log(0.5)) if False else np.log(0.5); continue
        i=i[ok]; j=j[ok]
        ea=np.exp(le[i]); eb=np.exp(le[j])
        if name=='circle': F=1/np.sqrt(1/ea**2+1/eb**2); F2=F
        else: F=H.s_star(ea,eb); F2=H.s_star(eb,ea)
        F=np.minimum(F,F2)
        if ovh:
            xmin=np.minimum(xs[i],xs[j])
            F=F*(1+2/xk+2.0**(-xmin))+2.0**(1-xmin)
        le[k]=min(np.log(0.5),np.min(np.log(F)))
    # report local slopes
    lx=np.log(xs); out=[]
    for X in report:
        kk=np.searchsorted(xs,X)
        if kk>=n: break
        k0=np.searchsorted(xs,X/10)
        sl=-(le[kk]-le[k0])/(lx[kk]-lx[k0])
        out.append(f"x=1e{np.log10(X):.0f}: e={np.exp(le[kk]):.3e} e*x={np.exp(le[kk])*X:.3g} locslope={sl:.3f}")
    return out
if __name__=='__main__':
    ovh=int(sys.argv[1]); 
    for name in sys.argv[2:]:
        print(name,'ovh',ovh); print('  '+'\n  '.join(run(name,ovh)),flush=True)
