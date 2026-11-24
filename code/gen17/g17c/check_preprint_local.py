# Sanity check: the preprint's local construction (Sec. 4.2-4.3: M0 = V0 f_I(W0), rows P_u(T) of (4.10))
# restricted to one swap-orbit with all Q pairs unequal, is an instance of the Lean setting
# (J(0)=X^{tensor l}, deg<=d, J(t)=K_t+E, K_t weakly t-nice), after relabelling columns by V0.
import itertools, sys
from fractions import Fraction as Fr
from math import comb

def rank(M):
    M=[r[:] for r in M]; rk=0; n=len(M); m=len(M[0]) if n else 0
    for c in range(m):
        p=None
        for i in range(rk,n):
            if M[i][c]!=0: p=i;break
        if p is None: continue
        M[rk],M[p]=M[p],M[rk]
        inv=1/M[rk][c]
        M[rk]=[v*inv for v in M[rk]]
        for i in range(n):
            if i!=rk and M[i][c]!=0:
                f=M[i][c]; M[i]=[a-f*b for a,b in zip(M[i],M[rk])]
        rk+=1
    return rk

def polymul(a,b):
    r=[Fr(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]+=x*y
    return r
def polyeval(a,t): return sum(c*Fr(t)**i for i,c in enumerate(a))

def run(Q,R,A0):
    l=Q; pts=list(itertools.product([0,1],repeat=l)); idx={x:i for i,x in enumerate(pts)}; n=len(pts)
    # f_I(z) = (-1)^A0 sum_{j<R} (-2)^j binom(z-A0, j)   (4.7); M0 = V0 f_I(W0), W0 = (1/2) sum (I - Pi_s)
    # Work in the character basis: eigenvalue of W0 on chi_z is |z|; M0 = H diag((-1)^|z| f(|z|)) H / 2^l
    def f(z):
        s=Fr(0)
        for j in range(R):
            b=Fr(1)
            for a in range(j): b*=Fr(z-A0-a)
            for a in range(1,j+1): b/=a
            s+=Fr(-2)**j*b
        return (-1)**A0*s
    lam=[(-1)**sum(z)*f(sum(z)) for z in pts]     # eigenvalues of M0
    def chi(z,x): return (-1)**sum(a*b for a,b in zip(z,x))
    M0=[[sum(lam[k]*chi(pts[k],x)*chi(pts[k],y) for k in range(n))/Fr(n) for y in pts] for x in pts]
    I=[[Fr(int(i==j)) for j in range(n)] for i in range(n)]
    # Lemma 4.2: M0[u,y]!=0 => |E_{u,y}| <= R-1, with E_{u,y} = {s : a_y(s) != d_u(s)} = {s : x_s == y_s}
    maxE=0
    for x in pts:
        for y in pts:
            if M0[idx[x]][idx[y]]!=0:
                maxE=max(maxE,sum(1 for s in range(l) if x[s]==y[s]))
    assert maxE<=R-1,(maxE,R)
    # P(T)[x,y] = (I-M0)[x,y] + M0[x,y] prod_{s in E}(1-T/s)   (4.10); points T = 1..Q (s is 1-based)
    def Pent(x,y):
        m=M0[idx[x]][idx[y]]
        p=[Fr(1)]
        for s in range(l):
            if x[s]==y[s]: p=polymul(p,[Fr(1),Fr(-1,s+1)])
        p=[m*c for c in p]
        p[0]+=I[idx[x]][idx[y]]-m
        return p
    comp=lambda y: tuple(1-b for b in y)
    # J = P * V0 (relabel columns y -> complement): J[x,y] = P[x, ~y]
    J={ (x,y):Pent(x,comp(y)) for x in pts for y in pts}
    deg=max(max([i for i,c in enumerate(p) if c!=0],default=0) for p in J.values())
    assert deg<=R-1
    # J(0) = X^{tensor l}
    for x in pts:
        for y in pts:
            assert polyeval(J[(x,y)],0)==(1 if y==comp(x) else 0)
    # common error E = (I - M0) V0
    E=[[I[idx[x]][idx[comp(y)]]-M0[idx[x]][idx[comp(y)]] for y in pts] for x in pts]
    rE=rank(E)
    assert rE==sum(1 for z in pts if not (A0<=sum(z)<=A0+R-1))
    tot01=0
    for t in range(1,Q+1):
        Jt=[[polyeval(J[(x,y)],t) for y in pts] for x in pts]
        K=[[Jt[i][j]-E[i][j] for j in range(n)] for i in range(n)]
        for x in pts:
            for y in pts:
                if x[t-1]!=y[t-1]: assert K[idx[x]][idx[y]]==0      # weak t-niceness of K_t
        B01=[[Jt[idx[x]][idx[y]] for y in pts if y[t-1]==1] for x in pts if x[t-1]==0]
        tot01+=rank(B01)
    d=deg; N=Q
    bound=sum(max(0,N-d-2*sum(x)) for x in pts)
    print(f"Q=l={Q} R={R} I=[{A0},{A0+R-1}] deg={deg} maxE={maxE} rank(E)={rE} (rho={rE}/{n}) "
          f"sum_k rank E = {Q*rE}; sum_k rank B01_k = {tot01}; Lean bound sum_x max(0,N-d-2|x|) = {bound}  "
          f"{'OK' if tot01>=bound and Q*rE>=bound else 'VIOLATION'}")

for (Q,R) in [(2,2),(3,2),(4,2),(4,3),(5,2),(5,3),(5,4),(6,3),(6,4)]:
    A0=(Q-R+1)//2
    run(Q,R,A0)
