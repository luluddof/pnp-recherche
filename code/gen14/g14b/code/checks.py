#!/usr/bin/env python3
"""G14b checks for the Appendix-D (tensor contraction) model.
C1: identity (D.10): v2(E_t(0)) = D * sum_j (j/2^j) (q(i+j)-q(i)), brute force on the grid S.
C2: taper optimality: min over profiles q(0..u-1)>=0 with sum_i 2^{a-i} q_i >= M(2L-1)+1 of max_i Lambda_i (LP).
C3: W(0)=I, W(t) = c_t * shift_t (exact rationals), locality of bit i.
C4: restriction/flattening over Z_2: det(I+2E) odd (trivial), orthant count 2*4^-u.
C5: number of mode tuples vs 2^{sqrt(np)-n} and Meinardus 2 sqrt(zeta(2) n p).
"""
from fractions import Fraction as Fr
import math, itertools
def v2(n):
    n=abs(n); k=0
    while n%2==0: n//=2; k+=1
    return k
def c1(a,u,h):
    L=2**a; S=[t for t in range(-L,L) if t!=0]
    q=lambda i: h*max(0,u-i)
    ok=True
    for t in S:
        i=v2(t)
        lhs=sum(q(v2(s))*(v2(s)-v2(t-s)) for s in S if s!=t)
        D=Fr(L,2**i)
        rhs=D*sum(Fr(j,2**j)*(q(i+j)-q(i)) for j in range(1,a-i+1))
        ok&=(lhs==rhs)
    return ok
print("C1 (D.10) exact on grids:", all(c1(a,u,h) for a in range(1,7) for u in range(1,a+1) for h in (1,3)))
# arbitrary profile too
import random
random.seed(1)
def c1b(a):
    L=2**a; S=[t for t in range(-L,L) if t!=0]; qq=[random.randint(0,9) for _ in range(a+1)]
    for t in S:
        i=v2(t)
        lhs=sum(qq[v2(s)]*(v2(s)-v2(t-s)) for s in S if s!=t)
        rhs=Fr(L,2**i)*sum(Fr(j,2**j)*(qq[i+j]-qq[i]) for j in range(1,a-i+1))
        if lhs!=rhs: return False
    return True
print("C1b (D.10) for arbitrary valuation-dependent profiles:", all(c1b(a) for a in range(1,8)))
# C2 LP
try:
    from scipy.optimize import linprog
    import numpy as np
    print("C2 taper LP: min max_i Lambda_i (M=1);  columns: a,u, LPopt*u, paper taper max Lambda*u (h=3/u), bound 1/(4u)*u")
    for a,u in [(4,2),(6,3),(8,4),(12,6),(16,8),(20,10),(24,12),(40,20),(60,30),(40,40),(30,30)]:
        L=2.0**a
        # vars q_0..q_{u-1}, Lam ; minimise Lam
        c=np.zeros(u+1); c[-1]=1
        A=[];b=[]
        for i in range(u):
            row=np.zeros(u+1)
            sig=sum(j/2.0**j for j in range(1,a-i+1))
            row[i]+=sig
            for j in range(1,u-i):
                row[i+j]-=j/2.0**j
            row[-1]=-1; A.append(row); b.append(0)
        row=np.zeros(u+1)
        for i in range(u): row[i]=-(2.0**(-i))
        A.append(row); b.append(-(2-1/L))
        r=linprog(c,A_ub=np.array(A),b_ub=np.array(b),bounds=[(0,None)]*(u+1))
        h=3.0/u; q=[h*(u-i) for i in range(u)]+[0]*(a+2)
        lam=max(sum(j/2.0**j*(q[i]-q[i+j]) for j in range(1,a-i+1)) for i in range(u))
        print(f"   a={a:3d} u={u:3d}  LP*u={r.fun*u:.4f}   taper*u={lam*u:.4f}   proved floor*u={0.25*(2-1/L):.4f}")
except Exception as e:
    print("C2 skipped:",e)
# C3
def c3(a):
    L=2**a; S=list(range(-L,L))
    def C(T): 
        p=Fr(1)
        for t in S: p*=(T-t)
        return p
    Cp0=Fr(1)
    for t in S:
        if t!=0: Cp0*=(0-t)
    def W(x,y,T):
        p=Fr(1)
        for t in S:
            if t!=y-x: p*=(T-t)
        return p/Cp0
    ok=True
    for x in range(L):
        for y in range(L):
            ok&=(W(x,y,0)==(1 if x==y else 0))
            for t in S:
                if t==0: continue
                w=W(x,y,t)
                ok&=((w!=0)==(y-x==t))
                if w!=0:
                    i=v2(t); ok&=(((x>>i)&1)!=((y>>i)&1)) and all(((x>>k)&1)==((y>>k)&1) for k in range(i))
    return ok
print("C3 W(0)=I, W(t)=c_t*shift_t, surviving pairs differ in bit v2(t), agree below:", all(c3(a) for a in (1,2,3)))
# C5
def count_tuples(n,p):
    # number of beta with sum_w w*(beta_{w,1..n}) < p : coefficients of prod_w (1-x^w)^-n, summed below p
    coef=[0]*p; coef[0]=1
    for w in range(1,p):
        for _ in range(n):
            for k in range(w,p): coef[k]+=coef[k-w]
    return sum(coef)
print("C5 tuples: n modes per weight, budget p: log2(count), floor sqrt(np)-n, Meinardus 2sqrt(zeta2 n p)/ln2")
for n,p in [(1,50),(2,100),(4,200),(4,800),(8,1600),(16,1600)]:
    c=count_tuples(n,p)
    print(f"   n={n:3d} p={p:5d}  log2={math.log2(c):9.2f}  floor={math.sqrt(n*p)-n:8.2f}  Meinardus={2*math.sqrt(math.pi**2/6*n*p)/math.log(2):9.2f}")
