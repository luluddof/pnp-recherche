# baseline constructions: tensor products and index deletion (restriction).  Writes frontier to cons_frontier.json (+ sols in cons/)
import numpy as np, json, os, itertools
from lib import *
os.makedirs('cons',exist_ok=True)
L=library()
by={}
def add(name,N,E,C):
    p,q=pq(N,E,C)
    by.setdefault(N,[]).append((p,q,name,E,C))
    return p,q
for name,N,E,C in L: add(name,N,E,C)
for N in range(1,13):
    for nm,E,C in trivial(N): add(nm,N,E,C)
# tensor products
for N in [4,6,8,9,10,12]:
    for N1 in range(2,N):
        if N%N1: continue
        N2=N//N1
        for (p1,q1,n1,E1,C1) in list(by.get(N1,[])):
            for (p2,q2,n2,E2,C2) in list(by.get(N2,[])):
                if q1*q2>2*N+8 or n1=='Swap' and n2=='Swap': continue
                if '(x)' in n1 and '(x)' in n2 and n1.count('(x)')>1: continue
                E=[np.kron(a,b) for a in E1 for b in E2]; C=[np.kron(a,b) for a in C1 for b in C2]
                p,q=pq(N,E,C)
                by.setdefault(N,[]).append((p,q,f"[{n1}]x[{n2}]",E,C))
# restriction: delete last index repeatedly (try every index) from N+1 down
for N in range(12,3,-1):
    for (p,q,n,E,C) in list(by.get(N+1,[])):
        bestp=None
        for i in range(N+1):
            keep=[j for j in range(N+1) if j!=i]
            E2=[e[np.ix_(keep,keep)] for e in E]; C2=[c[np.ix_(keep,keep)] for c in C]
            p2,q2=pq(N,E2,C2)
            if bestp is None or (p2,q2)<bestp[:2]: bestp=(p2,q2,E2,C2)
        by.setdefault(N,[]).append((bestp[0],bestp[1],f"del({n})",bestp[2],bestp[3]))
front={}
for N,lst in sorted(by.items()):
    b={}
    for p,q,n,E,C in lst:
        if q not in b or p<b[q][0]: b[q]=(p,n,E,C)
    # Pareto
    qs=sorted(b); out=[]; cur=10**9
    for q in qs:
        if b[q][0]<cur: cur=b[q][0]; out.append((q,b[q][0],b[q][1])); save(f"cons/c_{N}_{b[q][0]}_{q}.json",N,b[q][2],b[q][3],{"source":b[q][1]})
    front[N]=out
    print(N,[(p,q) for q,p,n in out if q<=2*N+8])
json.dump({str(k):v for k,v in front.items()},open('cons_frontier.json','w'),indent=0)
