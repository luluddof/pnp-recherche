import sys; sys.path.insert(0,'.')
from smith7 import *
# negative control: wrong alpha (2^a instead of 2^{a+1}) and wrong b-set (bit i differs)
a=4; L=16; bad=0; tot=0
for t in range(-L,L):
    if t==0 or v2(t)>=a: continue
    T=t+Fraction(2**a)
    C=Fraction(1)
    for q in range(-L,L): C*=(T-q)
    Cp0=Fraction(1)
    for q in range(-L,L):
        if q: Cp0*=-q
    i=v2(t); B=[x for x in range(L) if (x>>i)&1==0]
    E=[[C/(Cp0*(T-(y-x))) for y in B] for x in B]
    pred=predicted(a,i,1); N=max(pred)+30
    if any(e.denominator%2==0 for r in E for e in r): bad+=1; tot+=1; continue
    got=smith_vals([[red(e,N) for e in r] for r in E],N); tot+=1
    if got!=pred: bad+=1
print("wrong-alpha mismatches",bad,"of",tot)
