import numpy as np, itertools
k=4
P=['YIYX','IYZY','ZIII','IZYX','YYZI','XZXY','YYXX','XZYI','IXXZ','XIZY']
sg=[1,1,1,-1,1,-1,1,-1,1,-1]
def vec(p):
    a=b=0
    for i,ch in enumerate(p):
        if ch in 'XY': a|=1<<i
        if ch in 'ZY': b|=1<<i
    return a|(b<<k)
V=[vec(p) for p in P]
msk=(1<<k)-1
def Q(v): return bin((v&msk)&(v>>k)).count('1')&1
def B(u,v): return (bin((u&msk)&(v>>k)).count('1')+bin((v&msk)&(u>>k)).count('1'))&1
print('Q:',[Q(v) for v in V], 'sign*(-1)^Q:',[s*(-1)**Q(v) for s,v in zip(sg,V)])
G=np.array([[B(u,v) for v in V] for u in V]); print(G)
# rank of span
def rank(vs):
    vs=list(vs); r=0
    for bit in range(2*k):
        piv=next((i for i in range(r,len(vs)) if (vs[i]>>bit)&1),None)
        if piv is None: continue
        vs[r],vs[piv]=vs[piv],vs[r]
        for i in range(len(vs)):
            if i!=r and (vs[i]>>bit)&1: vs[i]^=vs[r]
        r+=1
    return r
print('rank span',rank(V))
# sums of triples that vanish / dependencies
deps=[c for r in range(3,6) for c in itertools.combinations(range(10),r) if np.bitwise_xor.reduce([V[i] for i in c])==0]
print('dependencies (size<=5):',deps)
