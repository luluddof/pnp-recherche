import sys,random,itertools; sys.path.insert(0,'.')
from smith7 import *
random.seed(3)
def det(M):
    M=[r[:] for r in M]; n=len(M); d=Fraction(1)
    for c in range(n):
        p=next((r for r in range(c,n) if M[r][c]!=0),None)
        if p is None: return Fraction(0)
        if p!=c: M[c],M[p]=M[p],M[c]; d=-d
        d*=M[c][c]
        for r in range(c+1,n):
            f=M[r][c]/M[c][c]
            if f: M[r]=[M[r][j]-f*M[c][j] for j in range(n)]
    return d
bad=0; tot=0
for a in (3,4,5):
    L=1<<a
    for t in [t for t in range(-L,L) if t and v2(t)<a]:
        for zeta in (Fraction(1),Fraction(6),Fraction(-1,3)):
            for b in (0,1):
                E,i=crossing(a,t,zeta,b); B=[x for x in range(L) if (x>>i)&1==b]
                T=t+Fraction(2**(a+1))*zeta
                C=Fraction(1)
                for q in range(-L,L): C*=(T-q)
                Cp0=Fraction(1)
                for q in range(-L,L):
                    if q: Cp0*=-q
                vc=v2f(C/Cp0); assert vc==a+1+v2f(zeta)
                for _ in range(6):
                    k=random.randint(1,len(B)); R=random.sample(range(len(B)),k); Cc=random.sample(range(len(B)),k)
                    d=det([[E[r][c] for c in Cc] for r in R])
                    xs=[B[r] for r in R]; ys=[B[c] for c in Cc]
                    form=k*vc+sum(v2(x-y) for x,y in itertools.combinations(xs,2))+sum(v2(x-y) for x,y in itertools.combinations(ys,2))-sum(min(i,v2(x-y)) for x in xs for y in ys)
                    tot+=1
                    if v2f(d)!=form: bad+=1; print("BAD",a,t,zeta,b,R,Cc)
print("minor-formula checks",tot,"bad",bad)
