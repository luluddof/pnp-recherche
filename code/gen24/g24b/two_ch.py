from objs import *
G=[[2,0],[0,2]]; xi=(Fr(1,2),Fr(1,2))
C2 = Copy(G, xi)
for (e,ha,hb) in [(0,0,0),(1,1,1),(2,1,1),(2,2,2),(1,2,2),(3,2,2),(2,3,3),(3,3,3),(4,2,2),(4,3,3),(4,4,4)]:
    R={g:rankQ(C2.B(g,ha,hb)) for g in C2.S}; L={g:rankQ(C2.O(g,e)) for g in C2.S}
    sc=max(min(R[g],L[g]) for g in C2.S)
    rk=flat_rank(C2,e,ha,hb)
    print("e=%d ha=%d hb=%d R=%s Lab=%s  maxSC=%d  sumSC=%d  loose=%d  sqrt(rkF)=%.2f rkF=%d"%(e,ha,hb,[R[g] for g in C2.S],[L[g] for g in C2.S],sc,sum(min(R[g],L[g]) for g in C2.S),min(len(C2.basis(ha)),len(C2.basis(hb)),len(C2.basis(e,True))),rk**.5,rk),flush=True)
