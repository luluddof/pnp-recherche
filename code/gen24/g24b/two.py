import sys, time
from objs import *
C1 = Copy([[2]], (Fr(1,2),))
C2 = Copy([[2,0],[0,2]], (Fr(1,2),Fr(1,2)))
print("S2 =", C2.S)
for (e,ha,hb) in [(0,0,0),(1,1,1),(2,0,0),(0,2,2),(2,1,1),(2,2,2),(1,2,2),(2,1,3),(3,1,1),(3,2,2),(2,3,3),(3,3,3),(4,2,2),(4,1,1)]:
    t=time.time()
    r1=flat_rank(C1,e,ha,hb); r2=flat_rank(C2,e,ha,hb); r2b=flat_rank(C2,e,ha,hb,2)
    dA=len(C2.basis(ha)); dM=len(C2.basis(e,True)); dB=len(C2.basis(hb))
    print("e=%d ha=%d hb=%d | 1copy rkF=%d | 2copy dims A,B,M=%d,%d,%d rkF=%d rkF2=%d sqrt=%.2f  (%.0fs)"%(e,ha,hb,r1,dA,dB,dM,r2,r2b,min(r2,r2b)**.5,time.time()-t),flush=True)
