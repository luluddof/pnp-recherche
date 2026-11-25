# Own LP + own rigorous verification (adaptive Bernstein enclosure in exact rationals) of Theorem-A certificates.
import numpy as np, sys
from fractions import Fraction as F
from math import comb
from scipy.optimize import linprog
def lp(q,a,k,marg=float(__import__("os").environ.get("MARG","2e-3")),grid=4000):
    js=list(range(1,k+1,2)); J=len(js)
    ys=np.unique(np.concatenate([np.linspace(0,1,grid),np.geomspace(1e-6,1,grid),np.linspace(1/q,1,grid)]))
    V=ys[:,None]**np.array(js)[None,:]; I=np.eye(J); sel=ys>=1/q
    A=np.vstack([np.hstack([I,-I]),np.hstack([-I,-I]),np.hstack([V,0*V]),np.hstack([-V,0*V]),np.hstack([-V[sel],0*V[sel]])])
    b=np.concatenate([np.zeros(2*J),(1-marg)*np.ones(len(ys)),(1-marg)*np.ones(len(ys)),-(a+marg)*np.ones(sel.sum())])
    r=linprog(np.r_[np.zeros(J),np.ones(J)],A_ub=A,b_ub=b,bounds=[(None,None)]*J+[(0,None)]*J,method='highs')
    assert r.status==0,r.message
    return js,r.x[:J]
def poly_coeffs(js,g,den=10**8):
    k=max(js); c=[F(0)]*(k+1)
    for j,v in zip(js,g): c[j]=F(round(v*den),den)
    return c
def shift_scale(c,lo,w):  # coefficients of p(lo + w t) in t
    n=len(c)-1; d=list(c)
    for i in range(n):          # Taylor shift by lo (Horner-like)
        for j in range(n-1,i-1,-1): d[j]+=lo*d[j+1]
    return [d[i]*w**i for i in range(n+1)]
def bern(c):  # power basis on [0,1] -> Bernstein
    n=len(c)-1
    return [sum(F(comb(i,j),comb(n,j))*c[j] for j in range(i+1)) for i in range(n+1)]
def min_on(c,lo,hi,depth=0):   # rigorous lower bound check machinery: returns min of Bernstein coeffs recursively until >= target
    pass
def check(c,lo,hi,pred,maxdepth=40):
    """verify pred(lower,upper) for enclosure of p on [lo,hi] by adaptive bisection; pred returns True(ok)/False(fail)/None(refine)"""
    stack=[(F(lo),F(hi),0)]; nint=0
    while stack:
        l,h,dp=stack.pop(); b=bern(shift_scale(c,l,h-l)); lb,ub=min(b),max(b); nint+=1
        r=pred(lb,ub)
        if r is True: continue
        # exact endpoint values are attained: if endpoint violates -> real failure
        if r is False or dp>maxdepth: return False,nint
        m=(l+h)/2; stack+= [(l,m,dp+1),(m,h,dp+1)]
    return True,nint
if __name__=='__main__':
    q=int(sys.argv[1]); a=F(sys.argv[2]); k=int(sys.argv[3])
    js,g=lp(q,float(a),k); c=poly_coeffs(js,g)
    L=sum(abs(x) for x in c)
    ok1=check(c,0,1,lambda lb,ub: True if (lb>=-1 and ub<=1) else None)
    ok2=check(c,F(1,q),1,lambda lb,ub: True if lb>=a else None)
    import sympy as sp
    y=sp.symbols('y'); P=sp.Poly(sum(sp.Rational(x.numerator,x.denominator)*y**i for i,x in enumerate(c)),y)
    s1=P.sub(sp.Poly(1,y)).count_roots(0,1)==0 and (P+1).count_roots(0,1)==0 and (P-sp.Rational(a.numerator,a.denominator)).count_roots(sp.Rational(1,q),1)==0
    print(f'q={q} a={float(a)} deg={k}: |g|<=1 on [0,1]: {ok1}  g>=a on [1/q,1]: {ok2}  sympy-sturm:{s1}  L={float(L):.4f}  => p >= {float(a/(1+a)):.4f} N^2 - {float(L/(1+a)):.3f} N')
    import json; json.dump({'q':q,'a':str(a),'coeffs':[str(x) for x in c],'L':str(L),'ok':bool(ok1[0] and ok2[0] and s1)},open(f'cert_q{q}_a{float(a)}_d{k}.json','w'))
