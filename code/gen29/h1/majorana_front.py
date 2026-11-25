# Majorana family: A = 2k+1 pairwise-anticommuting vectors (Jordan-Wigner) in F2^{2k}; y_s(w) = Q(s)+B(s,w) mod 2;
# J(w) = sum_s y_s (odd); CHECK beta(w) = F(J mod 4).  phi = f(u), u=(J-1)/2, f = Lagrange interpolant of F on a window of
# consecutive u values.  p exact (integer u), q = nnz of WHT(phi) (float, tol 1e-9 relative to rational structure).
import numpy as np, sys
from fractions import Fraction
from bell import popc, fwht
k=int(sys.argv[1]); n=4**k; msk=(1<<k)-1
A=[]
for j in range(k):
    for P in 'XY':
        a=1<<j; b=((1<<j)-1) | ((1<<j) if P=='Y' else 0); A.append(a|(b<<k))
t=0
for v in A: t^=v
A.append(t)
w=np.arange(n)
def Qv(x): return popc((x&msk)&(x>>k))&1
def Bv(x,s): return (popc((x&msk)&(s>>k)) + popc((s&msk)&(x>>k)))&1
J=np.zeros(n,np.int64)
for s in A: J+=(Qv(np.array([s]))[0] + Bv(w,s))&1
beta=1-2*Qv(w)
par=int(J[0]%2); assert np.all(J%2==par)
F={}
for j4 in (par,par+2):
    vals=set(beta[J%4==j4]); assert len(vals)==1, vals; F[j4]=vals.pop()
u=(J-par)//2
Fu=lambda uu: F[(2*uu+par)%4]
pts=[]
for lo in range(k+1):
    for hi in range(lo,k+1):
        W=list(range(lo,hi+1))
        def f(x):
            s=Fraction(0)
            for i in W:
                tt=Fraction(Fu(i))
                for jj in W:
                    if jj!=i: tt*=Fraction(x-jj,i-jj)
                s+=tt
            return s
        fv=[f(x) for x in range(k+1)]
        phi=np.array([float(v) for v in fv])[u]
        p=int(sum(np.sum(u==x) for x in range(k+1) if fv[x]!=Fu(x)))
        q=int(np.sum(np.abs(fwht(phi))/n>1e-9))
        pts.append((q,p,lo,hi))
pts.sort(); front=[]; bp=10**18
for q,p,lo,hi in pts:
    if p<bp: front.append((q,p,lo,hi)); bp=p
print(f'k={k} F(J mod 4)={F}  Majorana frontier (q,p,window in u):', front, flush=True)
