# Negative control + random Cauchy-like tests for Theorem U: E_xy=c/(X_x-Y_y), X,Y random 2-adic integers.
# Hypothesis H: v2(c) >= v2(X_x-X_x') for all x!=x' (gamma in O). Compare Smith(M) with union Smith(Per_k), over Q_2 (shifted).
import random, sys, itertools
from fractions import Fraction as Fr
from lib12 import *
random.seed(5)
NB=200; MOD=1<<NB
def sv(Mx,shift):
    return [s-shift for s in smith_vals([[toZ(z*2**shift,MOD) for z in r] for r in Mx],NB)]
stats={True:[0,0],False:[0,0]}
for trial in range(120):
    n=random.choice([2,3,4,5])
    while True:
        X=[random.randint(0,63) for _ in range(n)]
        if len(set(X))==n: break
    Y=[random.randint(0,63)*2+1+ (0) for _ in range(n)]  # make X-Y nonzero: X arbitrary, Y odd... ensure nonzero
    if any(x==y for x in X for y in Y): continue
    vc=random.randint(0,8); c=Fr(2**vc*random.choice([1,3,5,-7]))
    H=all(vc>=v2(X[x]-X[xp]) for x in range(n) for xp in range(n) if x!=xp)
    E=[[c/(X[x]-Y[y]) for y in range(n)] for x in range(n)]
    sh=8*n
    sm=sv(Mmat(E),sh)
    un=[0]
    for k in range(1,n+1): un+=sv(compound(E,k,per),sh)
    un.sort()
    eq=(sm==un)
    stats[H][0 if eq else 1]+=1
    if H and not eq: print("VIOLATION under H",X,Y,vc,sm,un)
print("H true: [equal, differ] =",stats[True]," H false:",stats[False])
