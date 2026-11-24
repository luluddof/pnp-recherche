from objs import *
import math
for n in range(1,6):
    G=[[2 if i==j else 0 for j in range(n)] for i in range(n)]
    C=Copy(G,(Fr(1,2),)*n,coset=[Fr(1,2)]*n)
    for E0s in ([1,2] if n<=3 else [1]):
        e=Fr(E0s)-Fr(n,4)
        r1=flat_rank(C,e,0,0,1); r2=flat_rank(C,e,0,0,2)
        dM=len(C.basis(e,True)); Dle=sum(len(C.basis(Fr(k)-Fr(n,4),True)) for k in range(E0s+1)); Dlt=Dle-dM
        L=max(rankQ(C.O(g,e)) for g in C.S); R=[rankQ(C.B(g,0,0)) for g in C.S]
        Nlemma=max(N for N in range(1,200) if N*(N-1)//2<=Dle) if Dlt==1 else None
        print("n=%d E0s=%d H=0: dimM_E0=%d dimM_<=E0=%d dimM_<E0=%d | flattening N<=floor(sqrt(min rkF))=%d (rkF=%d,%d) | lemma N<=%s"%(n,E0s,dM,Dle,Dlt,int(math.isqrt(min(r1,r2))),r1,r2,Nlemma),flush=True)
