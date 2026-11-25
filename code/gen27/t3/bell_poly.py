# Bell-model "polynomial method" family: phi(x,y) = f(x.y), f = interpolant of (-1)^m on a window A of d+1 integers.
# exact p, q via direct WHT on F_2^{2k} (small k) and via closed formulas (any k).
import numpy as np, itertools, sys
from fractions import Fraction
from math import comb, log2
def interp(A):
    # returns f as function m -> Fraction (Lagrange)
    def f(m):
        s=Fraction(0)
        for i,a in enumerate(A):
            t=Fraction((-1)**a)
            for j,b in enumerate(A):
                if j!=i: t*=Fraction(m-b,a-b)
            s+=t
        return s
    return f
def formula(k,A):
    f=interp(A); d=len(A)-1
    vals=[f(m) for m in range(k+1)]
    p=sum(comb(k,m)*3**(k-m) for m in range(k+1) if vals[m]!=(-1)**m)
    # binomial coefficients c_j: f(m) = sum_j c_j C(m,j)  (Newton forward differences at 0)
    c=[]; diff=vals[:]
    for j in range(k+1):
        c.append(diff[0]); diff=[diff[i+1]-diff[i] for i in range(len(diff)-1)]
    q=0
    for u in range(k+1):
        F=sum(c[j]*Fraction(1,4**j)*comb(k-u,j-u) for j in range(u,k+1))
        if F!=0: q+=comb(k,u)*3**u
    return p,q
def direct(k,A):
    f=interp(A); n=4**k
    xs=np.array(list(itertools.product([0,1],repeat=k)))
    M=xs@xs.T  # m = x.y
    vals={m:f(m) for m in range(k+1)}
    # use exact scaled integers: multiply by common denominator
    from math import lcm
    D=1
    for v in vals.values(): D=lcm(D,v.denominator)
    phi=np.vectorize(lambda m:int(vals[m]*D))(M).astype(np.int64)
    beta=((-1)**M)*D
    p=int((phi!=beta).sum())
    # WHT over F_2^{2k}
    h=phi.reshape(-1).copy()
    hsz=1
    while hsz<n:
        h=h.reshape(-1,2,hsz); a=h[:,0,:]+h[:,1,:]; b=h[:,0,:]-h[:,1,:]; h=np.stack([a,b],1).reshape(-1); hsz*=2
    q=int((h!=0).sum())
    return p,q
if __name__=="__main__":
    for k in range(2,8):
        for d in range(0,k+1):
            # best window of d+1 consecutive integers by mass
            best=None
            for a0 in range(0,k-d+1):
                A=list(range(a0,a0+d+1))
                mass=sum(comb(k,m)*3**(k-m) for m in A)
                if best is None or mass>best[0]: best=(mass,A)
            A=best[1]
            pf,qf=formula(k,A); pd,qd=direct(k,A)
            N=2**k
            print(k,d,A[0],"formula",pf,qf,"direct",pd,qd,"OK" if (pf,qf)==(pd,qd) else "MISMATCH", " theta=%.3f T=%.3f"%(log2(max(pf,1))/k,log2(qf)/k))
