import time
from objs import *
from pred2 import bound
C3 = Copy([[2,0,0],[0,2,0],[0,0,2]], (Fr(1,2),)*3)
print("S3 size", len(C3.S), flush=True)
for (e,ha,hb) in [(0,0,0),(1,1,1),(2,1,1),(1,2,2),(2,2,2),(3,1,1),(3,2,2)]:
    t=time.time()
    R={g:rankQ(C3.B(g,ha,hb)) for g in C3.S}; L={g:rankQ(C3.O(g,e)) for g in C3.S}
    sc=max(min(R[g],L[g]) for g in C3.S)
    rk=flat_rank(C3,e,ha,hb)
    print("3copy e=%d ha=%d hb=%d dims(V_ha,M_e)=(%d,%d) maxSC=%d sqrt(rkF)=%.2f rkF=%d blockbound=%d (%.0fs)"%(e,ha,hb,len(C3.basis(ha)),len(C3.basis(e,True)),sc,rk**.5,rk,bound(3,e,ha,hb),time.time()-t),flush=True)
