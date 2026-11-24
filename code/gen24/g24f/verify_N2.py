# exact verification: n=2 copies, E0s=1, H=0, N=2: a_x=e^{alpha_x}, b_1=e^{alpha_2}, b_2=e^{alpha_1}, y_x=|-_x>, u_x=<+_{other}|  (signs fixed by search)
from build import setup
from objs import Fr
import itertools, numpy as np
C,A,Bb,Mb,Bs,Os=setup(2,1,0,0)
ia={p:i for i,(p,m) in enumerate(A)}
im={(p,m):i for i,(p,m) in enumerate(Mb)}
def vec(d,idx):
    v=[Fr(0)]*d; v[idx]=Fr(1); return v
e1=(1,0); e2=(0,1)
a=[vec(len(A),ia[tuple(Fr(x) for x in e1)]),vec(len(A),ia[tuple(Fr(x) for x in e2)])]
b=[a[1],a[0]]
minus=[im[((Fr(-1,2),Fr(1,2)),())],im[((Fr(1,2),Fr(-1,2)),())]]
plus=[im[((Fr(3,2),Fr(1,2)),())],im[((Fr(1,2),Fr(3,2)),())]]
for s1,s2 in itertools.product([1,-1],repeat=2):
    Y=[vec(len(Mb),minus[0]),vec(len(Mb),minus[1])]
    U=[[Fr(s1)*z for z in vec(len(Mb),plus[1])],[Fr(s2)*z for z in vec(len(Mb),plus[0])]]
    ok=True
    for x,v,u,y in itertools.product(range(2),repeat=4):
        tot=Fr(0)
        for g in C.S:
            beta=sum(a[x][i]*Bs[g][i][j]*b[v][j] for i in range(len(A)) for j in range(len(Bb)))
            om=sum(U[u][k]*Os[g][k][l]*Y[y][l] for k in range(len(Mb)) for l in range(len(Mb)))
            tot+=beta*om
        if tot!=(1 if (x==v and u==y) else 0): ok=False
    print("signs",s1,s2,"identity exact:",ok)
