from objs import *
C = Copy([[2]], (Fr(1,2),))
print("S =", C.S)
for e in range(0,5):
  for H in range(0,9,2):
    for ha in [H//2]:
      hb=H-ha
      Bs={g:C.B(g,ha,hb) for g in C.S}; Os={g:C.O(g,e) for g in C.S}
      F=None; F2=None
      for g in C.S:
        t=kron(Bs[g],Os[g]); t2=kron(Bs[g],transpose(Os[g]))
        F=t if F is None else madd(F,t); F2=t2 if F2 is None else madd(F2,t2)
      R={g:rankQ(Bs[g]) for g in C.S}; L={g:rankQ(Os[g]) for g in C.S}
      dA=len(C.basis(ha)); dB=len(C.basis(hb)); dM=len(C.basis(e,True))
      rF=rankQ(F); rF2=rankQ(F2)
      print("e=%d ha=%d hb=%d dims(A,B,M)=(%d,%d,%d)  R=%s Lab=%s  sum R*Lab=%d  rkF=%d rkF2=%d  sqrt=%.2f  minRL=%s"%(e,ha,hb,dA,dB,dM,[R[g] for g in C.S],[L[g] for g in C.S],sum(R[g]*L[g] for g in C.S),rF,rF2,min(rF,rF2)**.5,[min(R[g],L[g]) for g in C.S]),flush=True)
