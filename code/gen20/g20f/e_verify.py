# Independent exact verification of crossing decompositions X = R + C (g19b json files).
import json, sys, itertools
from fractions import Fraction
def rank_exact(M):
    M=[[Fraction(v) for v in row] for row in M]; r=0; rows=len(M); cols=len(M[0]) if rows else 0
    for c in range(cols):
        piv=next((i for i in range(r,rows) if M[i][c]!=0),None)
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]
        for i in range(rows):
            if i!=r and M[i][c]!=0:
                f=M[i][c]/M[r][c]; M[i]=[a-f*b for a,b in zip(M[i],M[r])]
        r+=1
    return r
def rank_mod(M,p):
    M=[[int(v)%p for v in row] for row in M]; r=0; rows=len(M); cols=len(M[0])
    for c in range(cols):
        piv=next((i for i in range(r,rows) if M[i][c]),None)
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]; inv=pow(M[r][c],p-2,p); M[r]=[a*inv%p for a in M[r]]
        for i in range(rows):
            if i!=r and M[i][c]:
                f=M[i][c]; M[i]=[(a-f*b)%p for a,b in zip(M[i],M[r])]
        r+=1
    return r
for fn in sys.argv[1:]:
    d=json.load(open(fn)); N=d['N']; G=[[[int(Fraction(v)) for v in row] for row in g] for g in d['G']]; Gp=[[[int(Fraction(v)) for v in row] for row in g] for g in d['Gp']]
    q=len(G); assert len(Gp)==q
    R=range(N)
    C={}; X={}
    for x,y,u,v in itertools.product(R,R,R,R):
        C[x,y,u,v]=sum(G[l][x][u]*Gp[l][y][v] for l in range(q)); X[x,y,u,v]=int(x==v and y==u)
    Rt={k:X[k]-C[k] for k in X}
    rr=[[Rt[x,y,u,v] for u in R for v in R] for x in R for y in R]       # (x,y)|(u,v)
    rcC=[[C[x,y,u,v] for y in R for v in R] for x in R for u in R]       # (x,u)|(y,v)
    # also: maximum entry, entries in {0,+-1}
    vals=set(v for g in G+Gp for row in g for v in row)
    out=[fn.split('/')[-1],'N',N,'claimed',(d['p'],d['q']),'rr(R)=',rank_exact(rr),'rc(C)=',rank_exact(rcC),'entries',sorted(vals)]
    for p in (2,3,5,7):
        out+=[f'mod{p}:',(rank_mod(rr,p),rank_mod(rcC,p))]
    print(*out)
