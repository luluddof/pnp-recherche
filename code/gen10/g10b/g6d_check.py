# Independent re-check of the G6d determinant lemma (own code).
import itertools, random
from fractions import Fraction as Fr
from math import comb

def rank(M):
    M=[row[:] for row in M]; r=0
    if not M: return 0
    nc=len(M[0])
    for c in range(nc):
        p=None
        for i in range(r,len(M)):
            if M[i][c]!=0: p=i;break
        if p is None: continue
        M[r],M[p]=M[p],M[r]
        inv=1/M[r][c]
        for i in range(r+1,len(M)):
            if M[i][c]!=0:
                f=M[i][c]*inv
                M[i]=[a-f*b for a,b in zip(M[i],M[r])]
        r+=1
        if r==len(M): break
    return r

def bits(x,l): return [(x>>j)&1 for j in range(l)]

def minrankE(Jt,l,k):
    # min rank of E = J - K over K with K[x,y]=0 when x_k != y_k  = max(rank B, rank C)
    n=1<<l
    r0=[x for x in range(n) if not (x>>k)&1]; r1=[x for x in range(n) if (x>>k)&1]
    B=[[Jt[x][y] for y in r1] for x in r0]; C=[[Jt[x][y] for y in r0] for x in r1]
    return max(rank(B),rank(C))

def check_lemma(Jfun,l,d,points):
    """Jfun(t)-> 2^l x 2^l matrix; points = list of (tau, direction). Checks, for all S,
    (i) sum_k nullity(M(tau_k)) <= d 2^(l-s) on block A'xB', (ii) sum minrank E_k >= (N_S-d)2^(l-s)."""
    n=1<<l
    Jv={t:Jfun(t) for t,_ in points}
    mr={(t,c):minrankE(Jv[t],l,c) for t,c in points}
    worst=None
    for s in range(1,l+1):
        for S in itertools.combinations(range(l),s):
            mask=sum(1<<j for j in S)
            A=[x for x in range(n) if x&mask==0]; B=[y for y in range(n) if y&mask==mask]
            pts=[(t,c) for t,c in points if c in S]
            nul=0; lhs=0
            for t,c in pts:
                M=[[Jv[t][x][y] for y in B] for x in A]
                rk=rank(M); nul+=len(A)-rk
                assert mr[(t,c)]>=rk
                lhs+=mr[(t,c)]
            assert nul<=d*len(A),(S,nul)
            rhs=(len(pts)-d)*len(A)
            assert lhs>=rhs,(S,lhs,rhs)
            sl=lhs-rhs
            if rhs>0 and (worst is None or sl<worst[0]): worst=(sl,S,lhs,rhs)
    return mr,worst

def hadamard_conj(diag,l):
    # J = H diag H, H = 2^{-l/2} (-1)^{x.z}
    n=1<<l
    return [[sum(diag[z]*(-1)**(bin(x&z).count('1')+bin(y&z).count('1')) for z in range(n))/n for y in range(n)] for x in range(n)]

def chi(z): return (-1)**bin(z).count('1')

# sanity: H diag(chi) H = X^{(x)l}
for l in range(1,5):
    n=1<<l
    J0=hadamard_conj([Fr(chi(z)) for z in range(n)],l)
    assert all(J0[x][y]==(1 if y==(n-1)^x else 0) for x in range(n) for y in range(n))
print("H diag(chi) H = X^{(x)l}: ok")

# (1) G6d explicit d=1 family
def d1family(l,tau):
    n=1<<l
    def r(w):
        b=bits(w,l)
        return 1/tau[0] + (-1)**b[0]*(-1)**sum(b[1:])*sum(b[j]*(2/tau[j]-2/tau[0]) for j in range(1,l))
    def g(z,t): return chi(z)-t*chi(z)*r(z)
    return g
for l in range(2,6):
    tau=[Fr(k+1) for k in range(l)]
    g=d1family(l,tau); n=1<<l
    Jf=lambda t: hadamard_conj([g(z,t) for z in range(n)],l)
    pts=[(tau[k],k) for k in range(l)]
    mr,worst=check_lemma(Jf,l,1,pts)
    bad=[sum(1 for z in range(n) if not (z>>k)&1 and g(z,tau[k])!=g(z^(1<<k),tau[k])) for k in range(l)]
    print("d=1 family l=%d: minrank E_k=%s  #bad edges=%s  tightest (slack,S,lhs,rhs)=%s"%(l,[mr[p] for p in pts],bad,worst))
    assert [mr[p] for p in pts]==bad

# (2) G7a piecewise family
def lagrange(xs,ys):
    def p(s):
        tot=Fr(0)
        for i,(xi,yi) in enumerate(zip(xs,ys)):
            term=Fr(yi)
            for j,xj in enumerate(xs):
                if j!=i: term*=Fr(s-xj)/(xi-xj)
            tot+=term
        return tot
    return p
def piecewise(l,d,tau):
    blocks=[list(range(a,min(a+d+1,l+1))) for a in range(0,l+1,d+1)]
    pB={}
    for B in blocks:
        p=lagrange([l-2*w for w in B],[(-1)**w for w in B])
        for w in B: pB[w]=p
    def g(z,t):
        b=bits(z,l); s=[(-1)**x for x in b]; w=sum(b)
        return chi(z)+pB[w](sum((1-t/tau[j])*s[j] for j in range(l)))-pB[w](sum(s))
    return g
for l in range(2,6):
    for d in range(1,l):
        tau=[Fr(2*k+3,7) for k in range(l)]
        g=piecewise(l,d,tau); n=1<<l
        assert all(g(z,0)==chi(z) for z in range(n))
        Jf=lambda t: hadamard_conj([g(z,t) for z in range(n)],l)
        pts=[(tau[k],k) for k in range(l)]
        mr,worst=check_lemma(Jf,l,d,pts)
        bad=[sum(1 for z in range(n) if not (z>>k)&1 and g(z,tau[k])!=g(z^(1<<k),tau[k])) for k in range(l)]
        assert [mr[p] for p in pts]==bad
        print("piecewise l=%d d=%d: rank E_k=%s tightest=%s"%(l,d,bad,worst))

# (3) non-commutative monomial class, several points per direction:
# J(t) = X^{(x)l} * diag_y( prod_{tau in R_y} (1 - t/tau) ), |R_y| <= d  (adversarially many zero columns)
random.seed(1)
for trial in range(300):
    l=random.randint(1,4); d=random.randint(0,3); n=1<<l
    N=random.randint(1,2*l+2)
    taus=random.sample(range(1,40),N); points=[(Fr(t),random.randrange(l)) for t in taus]
    R=[random.sample(taus,min(d,N)) for y in range(n)]
    L=[[Fr(random.randint(-2,2)) if random.random()<0.3 else Fr(0) for _ in range(n)] for _ in range(n)]
    def Jf(t):
        col=[1]*n
        for y in range(n):
            v=Fr(1)
            for r in R[y]: v*=(1-Fr(t)/r)
            col[y]=v
        # plus a random degree<=d perturbation vanishing at 0 that is 'nice' nowhere in particular
        return [[(col[y] if y==(n-1)^x else 0) + (t**d)*L[x][y]*(1 if (trial%2 and d>=1) else 0) for y in range(n)] for x in range(n)]
    check_lemma(Jf,l,d,points)
print("monomial / perturbed non-commutative tests: 300 ok")

# (4) nullity-vs-order sharpness: M(t) = P diag(prod (1-t/tau)^..) Q with constant invertible P,Q
for trial in range(100):
    m=random.randint(1,5); d=random.randint(1,3)
    taus=random.sample(range(1,30),random.randint(1,6))
    Rj=[[random.choice(taus) for _ in range(d)] for _ in range(m)]
    while True:
        P=[[Fr(random.randint(-2,2)) for _ in range(m)] for _ in range(m)]
        Q=[[Fr(random.randint(-2,2)) for _ in range(m)] for _ in range(m)]
        if rank(P)==m and rank(Q)==m: break
    tot=0
    for t in taus:
        D=[Fr(1)]*m
        for j in range(m):
            for r in Rj[j]: D[j]*=(1-Fr(t)/r)
        M=[[sum(P[i][k]*D[k]*Q[k][j] for k in range(m)) for j in range(m)] for i in range(m)]
        tot+=m-rank(M)
    assert tot<=m*d
print("sum nullity <= m d (incl. repeated roots): ok")
