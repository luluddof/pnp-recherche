from trank import *
from fractions import Fraction as F
G=[[2]]; xi=(F(1,2),); g0=(0,); ga=(1,)
print("dims V'_E (xi=a/2):",[dimM(G,xi,E) for E in range(9)])
print("dims V_E conformal:",[dimM(G,(0,),E) for E in range(9)])
for H in range(2,11,2):
    h=H//2
    r0=R_channel(G,xi,g0,h,h)[0]; ra=R_channel(G,xi,ga,h,h)[0]
    r0s=[R_channel(G,xi,g0,ha,H-ha)[0] for ha in range(H+1)]
    ras=[R_channel(G,xi,ga,ha,H-ha)[0] for ha in range(H+1)]
    print("H=%d: R_0 splits %s ; R_alpha splits %s ; dimV'_{H/2}=%d"%(H,r0s,ras,dimM(G,xi,h)))
for E0 in range(0,8):
    l0=Lab_channel(G,xi,g0,E0)[0]; la,bl=Lab_channel(G,xi,ga,E0)
    print("E0=%d: dim M_E0=%d Lab_0=%d Lab_alpha=%d  conformal dim_E0=%d  blocks %s"%(E0,dimM(G,xi,E0),l0,la,dimM(G,(0,),E0),[(int(b[0][0]),b[3],b[4],b[5]) for b in bl]))
