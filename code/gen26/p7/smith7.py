# Own independent check of Theorem 1: exact 2-adic Smith form of crossing block of (D.7)
import sys
from fractions import Fraction
def v2(n):
    n=abs(n); 
    if n==0: return 10**9
    return (n & -n).bit_length()-1
def v2f(f): return v2(f.numerator)-v2(f.denominator)
def red(f,N):  # Fraction with odd denominator -> int mod 2^N
    M=1<<N
    assert f.denominator%2==1
    return (f.numerator*pow(f.denominator,-1,M))%M
def smith_vals(A,N):
    # A: list of lists of ints mod 2^N; returns sorted valuations (N means 'zero mod 2^N')
    M=1<<N; A=[r[:] for r in A]; n=len(A); m=len(A[0]); out=[]
    rows=list(range(n)); cols=list(range(m))
    while rows and cols:
        best=None
        for r in rows:
            for c in cols:
                x=A[r][c]
                if x:
                    v=v2(x)
                    if best is None or v<best[0]:
                        best=(v,r,c)
                        if v==0: break
            if best and best[0]==0: break
        if best is None: break
        v,r,c=best; out.append(v)
        u=A[r][c]>>v; ui=pow(u,-1,M)
        for r2 in rows:
            if r2==r: continue
            x=A[r2][c]
            if x:
                f=((x>>v)*ui)%M
                Ar=A[r]; A2=A[r2]
                for c2 in cols: A2[c2]=(A2[c2]-f*Ar[c2])%M
        rows.remove(r); cols.remove(c)
    out+= [N]*(min(n,m)-len(out))
    return sorted(out)
def crossing(a,t,zeta,b):
    L=1<<a; S=range(-L,L); T=t+Fraction(2**(a+1))*zeta
    C=Fraction(1)
    for q in S: C*=(T-q)
    Cp0=Fraction(1)
    for q in S:
        if q!=0: C0=Cp0; Cp0*=(0-q)
    i=v2(t)
    B=[x for x in range(L) if (x>>i)&1==b]
    E=[[C/(Cp0*(T-(y-x))) for y in B] for x in B]
    return E,i
def predicted(a,i,zeta):
    vz=v2f(Fraction(zeta)) if zeta!=0 else None
    out=[]
    for k in range(2**(a-i-1)):
        vk=sum(v2(j) for j in range(1,k+1))
        out+= [(a-i+1)+vz+2*k+2*vk]*(2**i)
    return sorted(out)
if __name__=="__main__":
    import random
    random.seed(7)
    amax=int(sys.argv[1]); nfail=0; ncase=0
    zetas=[Fraction(1),Fraction(-1),Fraction(2),Fraction(-3),Fraction(1,3),Fraction(-7,5),Fraction(8),Fraction(-12),Fraction(5,7)*4,Fraction(96)]
    for a in range(1,amax+1):
        L=1<<a
        ts=[t for t in range(-L,L) if t!=0 and v2(t)<a]
        if a>=6: ts=random.sample(ts,min(len(ts),12))+[-L+1,L-1,-(1<<(a-1)) if a>1 else -1, 3<<(a-2) if a>=2 else 1]
        for t in sorted(set(ts)):
            for zeta in (zetas if a<=4 else random.sample(zetas,3)):
                for b in (0,1):
                    E,i=crossing(a,t,zeta,b)
                    pred=predicted(a,i,zeta)
                    N=max(pred)+20
                    got=smith_vals([[red(e,N) for e in row] for row in E],N)
                    ncase+=1
                    if got!=pred:
                        nfail+=1; print("FAIL",a,t,zeta,b,got,pred)
        print("a",a,"cases so far",ncase,"fails",nfail,flush=True)
