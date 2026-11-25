# Exact certificate for Theorem A at fixed q: odd poly g with rational coeffs, verified exactly:
#   g(y) <= 1 on [0,1],  g(y) >= a on [1/q,1]   =>  p >= (a N^2 - L N)/(1+a),  L = sum |g_j|.
import numpy as np, sys, sympy as sp
from fractions import Fraction
from scipy.optimize import linprog
def solve(q,a,k,marg=1e-3,grid=3000):
    js=list(range(1,k+1,2)); J=len(js)
    ys=np.unique(np.concatenate([np.linspace(0,1,grid),np.geomspace(1e-5,1,grid)]))
    V=np.array([[y**j for j in js] for y in ys]); I=np.eye(J)
    c=np.concatenate([np.zeros(J),np.ones(J)]); sel=ys>=1.0/q
    A=np.vstack([np.hstack([I,-I]),np.hstack([-I,-I]),np.hstack([V,0*V]),np.hstack([-V,0*V]),np.hstack([-V[sel],0*V[sel]])])
    b=np.concatenate([np.zeros(2*J),(1-marg)*np.ones(len(ys)),(1-marg)*np.ones(len(ys)),-(a+marg)*np.ones(sel.sum())])
    r=linprog(c,A_ub=A,b_ub=b,bounds=[(None,None)]*J+[(0,None)]*J,method='highs')
    assert r.status==0, r.message
    return js,r.x[:J]
def verify(q,a,js,g,den=10**9):
    y=sp.symbols('y'); coeffs=[sp.Rational(int(round(v*den)),den) for v in g]
    P=sum(c*y**j for c,j in zip(coeffs,js)); a=sp.Rational(a)
    def maxmin(expr,lo,hi):
        # exact: extremes at endpoints or real roots of derivative in [lo,hi]
        d=sp.Poly(sp.diff(expr,y),y); pts=[lo,hi]
        for iv in d.intervals(): 
            l,h=iv[0]
            if h>=lo and l<=hi: pts.append(sp.Rational(max(l,lo)+min(h,hi),2))  # bracket
        return pts
    # rigorous: check g-1 has no root in (0,1] with sign change and g(1)<=1 ... use Sturm counts
    ok1 = sp.Poly(P-1,y).count_roots(0,1)==0 and P.subs(y,sp.Rational(1,2))<1
    ok3 = sp.Poly(P+1,y).count_roots(0,1)==0 and P.subs(y,sp.Rational(1,2))>-1
    ok2 = sp.Poly(P-a,y).count_roots(sp.Rational(1,q),1)==0 and P.subs(y,1)>a
    L=sum(abs(c) for c in coeffs)
    return ok1 and ok3,ok2,L
if __name__=='__main__':
    q=int(sys.argv[1]); a=float(sys.argv[2]); k=int(sys.argv[3])
    js,g=solve(q,a,k); ok1,ok2,L=verify(q,a,js,g)
    print('q',q,'a',a,'deg',k,'verified |g|<=1 on [0,1]:',ok1,' g>=a on [1/q,1]:',ok2,' L=%.4f'%float(L))
    print(' => p >= (%.4f N^2 - %.4f N)/%.4f = %.4f N^2 - %.4f N'%(a,float(L),1+a,a/(1+a),float(L)/(1+a)))
