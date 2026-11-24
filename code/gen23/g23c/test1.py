from trank import *
from fractions import Fraction as F
print([len(fock(1,k)) for k in range(8)], [len(fock(2,k)) for k in range(6)])
# App B: Heisenberg (charge 0 only)
G=[[1000]]; xi=(0,)
for h in range(1,7):
    R,_=R_channel(G,xi,(0,),h,h); L,_=Lab_channel(G,xi,(0,),h)
    # also all splits ha+hb=2h
    Rs=[R_channel(G,xi,(0,),ha,2*h-ha)[0] for ha in range(0,2*h+1)]
    print("Heis h=%d: R_0(h,h)=%d Lab_0(h)=%d  p(h)=%d   R_0 over splits ha=0..2h:"%(h,R,L,len(fock(1,h))),Rs)
# A1 regraded xi=alpha/2
G=[[2]]; xi=(F(1,2),)
print("S_xi A1, xi=a/2:", S_xi(G,xi))
