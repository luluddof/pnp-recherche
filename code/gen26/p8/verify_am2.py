# Verify the explicit construction proving C1'-right (R_p(M(E)) >= square-free count) at level i = a-2.
# For each class rho (residue mod 2^i), rows l0=l_{rho,0}, l1=l_{rho,1} (x = rho + b 2^i + 2^{i+1} x'), u = (t+alpha zeta)/2^i (odd):
#   F_{rho,{}} = 1, F_{rho,{0}} = l0 (weight w0), F_{rho,{1}} = u l0 - (u+2) l1 (weight w0+2), F_{rho,{0,1}} = l0 l1 + 2^{w0} l0 (weight 2w0+1)
# F_T = prod_rho F_{rho,T_rho}. Check: 2^{w(T)} | F_T on the whole cube, and g_T = F_T/2^{w(T)} mod 2 form a basis of Fun(F_2^n,F_2).
import sys, itertools
from fractions import Fraction
K=200; MOD=1<<K
def v2(n):
    n=abs(n); return K if n==0 else (n&-n).bit_length()-1
def v2q(q): q=Fraction(q); return v2(q.numerator)-v2(q.denominator)
def toZ(fr): fr=Fraction(fr); return (fr.numerator*pow(fr.denominator,-1,MOD))%MOD
def run(a,t,zeta,b=0):
    L=1<<a; i=v2(t); assert i==a-2
    T=t+(1<<(a+1))*Fraction(zeta); Cp0=1
    for q in range(-L,L):
        if q: Cp0*=-q
    CT=Fraction(1)
    for q in range(-L,L): CT*=(T-q)
    B=[x for x in range(L) if (x>>i)&1==b]; idx={x:k for k,x in enumerate(B)}; n=len(B)
    E=[[toZ(CT/(Cp0*(T-(y-x)))) for y in B] for x in B]
    w0=3+v2q(zeta); assert v2q(CT/Cp0)==w0+i
    u=toZ(T/2**i)
    classes=range(1<<i)
    rows={(r,xp):idx[r+b*(1<<i)+(xp<<(i+1))] for r in classes for xp in (0,1)}
    pts=range(1<<n)
    # l_x(X)
    lin=[[0]*(1<<n) for _ in range(n)]
    for x in range(n):
        for Xm in range(1,1<<n):
            low=Xm&-Xm; k=low.bit_length()-1
            lin[x][Xm]=(lin[x][Xm^low]+E[x][k])%MOD
    fac={}
    for r in classes:
        l0=lin[rows[(r,0)]]; l1=lin[rows[(r,1)]]
        fac[(r,0)]=([1]*(1<<n),0)
        fac[(r,1)]=(l0,w0)
        fac[(r,2)]=([(u*p-(u+2)*q)%MOD for p,q in zip(l0,l1)],w0+2)
        fac[(r,3)]=([(p*q+(p<<w0))%MOD for p,q in zip(l0,l1)],2*w0+1)
    # check each factor's divisibility separately (sufficient), and collect reductions
    red={}
    ok=True
    for key,(f,w) in fac.items():
        if min(v2(z) for z in f)<w: ok=False; print("  factor",key,"min v2",min(v2(z) for z in f),"<",w)
        red[key]=[ (z>>w)&1 for z in f]
    # basis check: products of reductions, as bit-vectors over points
    vecs=[]; weights=[]
    for Ts in itertools.product(range(4),repeat=len(classes)):
        v=[1]*(1<<n); w=0
        for r,s in zip(classes,Ts):
            rv=red[(r,s)]; v=[p&q for p,q in zip(v,rv)]; w+=fac[(r,s)][1]
        vecs.append(int(''.join(map(str,v)),2)); weights.append(w)
    # GF(2) rank
    basis={}
    for v in vecs:
        while v:
            h=v.bit_length()-1
            if h in basis: v^=basis[h]
            else: basis[h]=v; break
    rank=len(basis)
    print(f"a={a} t={t} zeta={zeta} b={b} w0={w0} n={n}: factors divisible={ok}, rank={rank}/{1<<n}",flush=True)
    return ok and rank==(1<<n), sorted(weights), E
if __name__=='__main__':
    allok=True
    for a in (2,3,4):
        L=1<<a; i=a-2
        ts=[t for t in range(-L,L) if t!=0 and v2(t)==i]
        for t in ts[:6]:
            for zeta in (1,-3,Fraction(1,3),2,-12,Fraction(5,7)):
                for b in (0,1):
                    ok,ws,E=run(a,t,zeta,b); allok&=ok
    print("ALL OK" if allok else "FAILURES")
