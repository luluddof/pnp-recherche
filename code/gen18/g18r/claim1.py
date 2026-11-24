# Brute-force test of Proposition 1: 2^c <= |C{0,1}^m| <= min over row bases prod(||row||_1+1); mod p: <= p^{rank_p}
import numpy as np, itertools, random
from fractions import Fraction
def rank_Q(rows):
    M=[[Fraction(x) for x in r] for r in rows]; rk=0
    for c in range(len(M[0]) if M else 0):
        p=next((i for i in range(rk,len(M)) if M[i][c]!=0),None)
        if p is None: continue
        M[rk],M[p]=M[p],M[rk]
        for i in range(len(M)):
            if i!=rk and M[i][c]!=0:
                f=M[i][c]/M[rk][c]; M[i]=[a-f*b for a,b in zip(M[i],M[rk])]
        rk+=1
    return rk
def rank_p(rows,p):
    M=[[x%p for x in r] for r in rows]; rk=0
    for c in range(len(M[0]) if M else 0):
        pv=next((i for i in range(rk,len(M)) if M[i][c]),None)
        if pv is None: continue
        M[rk],M[pv]=M[pv],M[rk]; iv=pow(M[rk][c],-1,p); M[rk]=[(a*iv)%p for a in M[rk]]
        for i in range(len(M)):
            if i!=rk and M[i][c]:
                f=M[i][c]; M[i]=[(a-f*b)%p for a,b in zip(M[i],M[rk])]
        rk+=1
    return rk
rng=random.Random(5); viol=0; tests=0; tightlo=0; tightup=0; maxgap=0
for trial in range(6000):
    k=rng.randint(1,4); m=rng.randint(1,7); ent=rng.choice([1,2,5])
    C=[[rng.randint(-ent,ent) if rng.random()<0.6 else 0 for _ in range(m)] for _ in range(k)]
    if rng.random()<0.3 and k>1: C[-1]=[a+b for a,b in zip(C[0],C[1%k])]   # force dependencies
    c=rank_Q(C)
    X=np.array(list(itertools.product((0,1),repeat=m))); Y=X@np.array(C).T
    cnt=len(set(map(tuple,Y)))
    ub=None
    for I in itertools.combinations(range(k),c):
        if rank_Q([C[i] for i in I])==c or c==0:
            val=1
            for i in I: val*=sum(abs(x) for x in C[i])+1
            ub=val if ub is None else min(ub,val)
    nzc=sum(1 for j in range(m) if any(C[i][j] for i in range(k)))
    ok=(2**c<=cnt<=ub) and cnt<=2**nzc
    for p in (2,3,5):
        cp=len(set(map(tuple,Y%p))); ok&=(cp<=p**rank_p(C,p)) and rank_p(C,p)<=c
    tests+=1; viol+=(not ok); tightlo+=(cnt==2**c); tightup+=(cnt==ub)
print(f"tests={tests} violations={viol} lower-tight={tightlo} upper-tight={tightup}")
print("example C=(1,-1):", len({a-b for a in (0,1) for b in (0,1)}), "values; rank 1; bounds 2..3")
