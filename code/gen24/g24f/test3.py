# rank of the combined form B_0 + kappa*B_alpha on V'_h x V'_h (A1, xi=alpha/2), exact
from trank import *
from fractions import Fraction as F
G=[[2]]; xi=(F(1,2),)
def combined(h,kappa):
    r=1; rows=[];
    A=[(b,m) for b in lattice_points(G,2*(h+4)) for k in [h-wt(G,xi,b)] if k>=0 for m in fock(1,int(k))]
    Bs=A
    M=[]
    for (beta,am) in A:
        row=[]
        for (delta,bm) in Bs:
            g=beta[0]+delta[0]
            v=F(0)
            if g==0: v=Bval(G,beta,am,delta,bm)
            elif g==1: v=kappa*Bval(G,beta,am,delta,bm)
            row.append(v)
        M.append(row)
    return len(A), rank_q(M)
for h in range(0,7):
    print(h, combined(h,F(1)), combined(h,F(3,7)), "dim V'_h=",dimM(G,xi,h))
