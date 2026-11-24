# build exact B_g (pairs) and O_g (coset labels) for n copies of A1, xi=alpha/2, labels in V_{alpha/2+L}
# E0s = SHIFTED label weight (bottom of M = 0).
from objs import *
import numpy as np
def setup(n, E0s, ha, hb):
    G=[[2 if i==j else 0 for j in range(n)] for i in range(n)]
    C=Copy(G,(Fr(1,2),)*n,coset=[Fr(1,2)]*n)
    e=Fr(E0s)-Fr(n,4)
    A=C.basis(ha); Bb=C.basis(hb); Mb=C.basis(e,True)
    Bs={g:np.array(C.B(g,ha,hb),dtype=object) for g in C.S}
    Os={g:np.array(C.O(g,e),dtype=object) for g in C.S}
    return C,A,Bb,Mb,Bs,Os
if __name__=="__main__":
    import sys
    n,E0s,ha,hb=map(int,sys.argv[1:5])
    C,A,Bb,Mb,Bs,Os=setup(n,E0s,ha,hb)
    print("S",C.S,"dims",len(A),len(Bb),len(Mb))
    for g in C.S:
        print(g,"rkB",rankQ(Bs[g].tolist()) if len(A) else 0,"rkO",rankQ(Os[g].tolist()) if len(Mb) else 0)
    print(Mb)
