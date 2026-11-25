# D=4 "doubly-even" weight-statistic family.  C = random self-orthogonal code of dim K=2k, length n (projective: distinct
# nonzero columns); c chosen with c.x = wt(x)/2 mod 2 on C => all weights of c+C are = wt(c) mod 4.
# q(w) = (wt(c+Gw)-wt(c))/4 mod 2 is a quadratic function; if its polar form |x&y|/2 is nondegenerate, (-1)^q is affinely
# equivalent to +-beta (IP_k).  phi = f(u), u = (wt - j0)/4, f interpolates (-1)^u on a window.  p exact; q = WHT count.
import numpy as np, sys
from fractions import Fraction
from bell import fwht
k=int(sys.argv[1]); n=int(sys.argv[2]); tries=int(sys.argv[3]); seed=int(sys.argv[4]) if len(sys.argv)>4 else 0
K=2*k; rng=np.random.default_rng(seed)
def rank2(M):
    M=M.copy()%2; r=0; rows,cols=M.shape
    for c_ in range(cols):
        piv=np.where(M[r:,c_]==1)[0]
        if len(piv)==0: continue
        p=r+piv[0]; M[[r,p]]=M[[p,r]]
        for i in range(rows):
            if i!=r and M[i,c_]: M[i]^=M[r]
        r+=1
        if r==rows: break
    return r
def selforth():
    rows=[]
    for _ in range(4000):
        if len(rows)==K: break
        v=rng.integers(0,2,n)
        if v.sum()%2: continue
        if any((v&r_).sum()%2 for r_ in rows): continue
        cand=np.array(rows+[v])
        if rank2(cand)<len(rows)+1: continue
        rows.append(v)
    if len(rows)<K: return None
    G=np.array(rows)   # K x n generator (rows = basis codewords)
    cols=[tuple(G[:,j]) for j in range(n)]
    if len(set(cols))<n or any(sum(c_)==0 for c_ in cols): return None
    return G
best=None
W=np.array([[ (w>>i)&1 for i in range(K)] for w in range(4**k)],dtype=np.int64)
for t in range(tries):
    G=selforth()
    if G is None: continue
    # polar form on basis: B_ij = |g_i & g_j|/2 mod 2
    B=np.array([[((G[i]&G[j]).sum()//2)%2 if i!=j else 0 for j in range(K)] for i in range(K)])
    if rank2(B)<K: continue
    # c: solve c.g_i = wt(g_i)/2 mod 2  (any solution); find via least: c in row space of G^T... brute: solve G c = r
    r=(G.sum(1)//2)%2
    A=np.concatenate([G,r[:,None]],1).copy()
    # gaussian elimination to find c
    Mx=A.copy(); piv=[]; rr=0
    for col in range(n):
        ps=np.where(Mx[rr:,col]==1)[0]
        if len(ps)==0: continue
        p=rr+ps[0]; Mx[[rr,p]]=Mx[[p,rr]]
        for i in range(K):
            if i!=rr and Mx[i,col]: Mx[i]^=Mx[rr]
        piv.append(col); rr+=1
        if rr==K: break
    c=np.zeros(n,np.int64)
    for i,col in enumerate(piv): c[col]=Mx[i,n]
    code=(W@G)%2 ^ c[None,:]
    wt=code.sum(1); j0=int(wt[0]); assert np.all((wt-j0)%4==0)
    u=(wt-j0)//4; sgn=(1-2*(u%2))
    vals=sorted(set(u.tolist())); cnt={v:int((u==v).sum()) for v in vals}
    # frontier over windows
    pts=[]
    for a in range(len(vals)):
        for b2 in range(a,len(vals)):
            Wn=vals[a:b2+1]
            # interpolate (-1)^v on window values Wn (they may be non-consecutive integers)
            def f(x):
                s=Fraction(0)
                for i in Wn:
                    tt=Fraction((-1)**(i%2))
                    for j in Wn:
                        if j!=i: tt*=Fraction(x-j,i-j)
                    s+=tt
                return s
            fv={v:f(v) for v in vals}
            phi=np.array([float(fv[x]) for x in u])
            p=sum(cnt[v] for v in vals if fv[v]!=(-1)**(v%2))
            q=int(np.sum(np.abs(fwht(phi))/4**k>1e-9))
            pts.append((q,p))
    pts.sort(); fr=[]; bp=10**18
    for q,p in pts:
        if p<bp: fr.append((q,p)); bp=p
    print(f'k={k} n={n} try={t}: u-values {cnt}  frontier {fr}',flush=True)
    break
