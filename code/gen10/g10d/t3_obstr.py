"""G10d test 3: obstructions. (a) Prop H1 vs exact minimal decision trees incl. irregular graphs and a more
permissive tree model (queries of any vertex allowed, only consistent leaves counted);
(b) Gottlieb rank of disjointness matrices; (c) Prop S incidence bound on random regular graphs, arbitrary clusterings."""
import math, random, itertools, sys
import networkx as nx, numpy as np
from functools import lru_cache
from fractions import Fraction
sys.setrecursionlimit(10**6)
def F(D,j,memo={}):
    if j<=0: return 1
    if (D,j) not in memo: memo[(D,j)]=F(D,j-1)+F(D,j-D-1)
    return memo[(D,j)]
def minleaves(G,k):
    n=G.number_of_nodes(); adj=[0]*n
    for u,v in G.edges(): adj[u]|=1<<v; adj[v]|=1<<u
    @lru_cache(maxsize=None)
    def L(S):
        if bin(S).count('1')<=k: return 1
        best=None; T=S
        while T:
            b=T&-T; v=b.bit_length()-1; T^=b
            c=L(S&~b)+L(S&~b&~adj[v])
            if best is None or c<best: best=c
        return best
    return L((1<<n)-1)
rng=random.Random(5); ok=True; tight=0; cnt=0
graphs=[]
for n,d in [(12,3),(14,3),(14,4),(15,4),(16,5),(16,3)]:
    graphs.append(("reg%d_%d"%(n,d), nx.random_regular_graph(d,n,seed=rng.randrange(10**6))))
for n,p in [(13,0.2),(14,0.35),(15,0.15),(14,0.6)]:
    graphs.append(("gnp%d_%.2f"%(n,p), nx.gnp_random_graph(n,p,seed=rng.randrange(10**6))))
graphs += [("star13",nx.star_graph(12)),("path14",nx.path_graph(14)),("cycle15",nx.cycle_graph(15)),("K7",nx.complete_graph(7)),
           ("K44",nx.complete_bipartite_graph(4,4)),("petersen",nx.petersen_graph()),("empty10",nx.empty_graph(10)),
           ("5K3",nx.disjoint_union_all([nx.complete_graph(3)]*5)), ("3K5",nx.disjoint_union_all([nx.complete_graph(5)]*3))]
for name,G in graphs:
    G=nx.convert_node_labels_to_integers(G); n=G.number_of_nodes(); D=max(dict(G.degree()).values()) if n else 0
    for k in range(0,n+1):
        ml=minleaves(G,k); f=F(D,n-k); cnt+=1
        if ml<f: ok=False; print("VIOLATION",name,k,ml,f)
        if ml==f: tight+=1
    print(name,"n=%d Delta=%d"%(n,D),"k=0:",minleaves(G,0),">=",F(D,n),"  k=n//2:",minleaves(G,n//2),">=",F(D,n-n//2))
print("H1: %d (graph,k) cases, all >= F: %s ; equality in %d"%(cnt,ok,tight))
# closed-form lower bound F(j) >= x^(j-D)
for D in [1,2,3,5,8]:
    lo,hi=1.0,2.0
    for _ in range(200):
        x=(lo+hi)/2
        if D*math.log(x)+math.log(x-1)>0: hi=x
        else: lo=x
    assert all(F(D,j)>=x**(j-D)*(1-1e-9) for j in range(0,400)), D
print("F(j) >= x^(j-Delta) verified for Delta in {1,2,3,5,8}, j<400")

# (b) Gottlieb
def rank_exact(M):
    M=[[Fraction(x) for x in r] for r in M]; r=0; rows=len(M); cols=len(M[0]) if M else 0
    for c in range(cols):
        p=None
        for i in range(r,rows):
            if M[i][c]!=0: p=i; break
        if p is None: continue
        M[r],M[p]=M[p],M[r]
        for i in range(r+1,rows):
            if M[i][c]!=0:
                f=M[i][c]/M[r][c]
                M[i]=[a-f*b for a,b in zip(M[i],M[r])]
        r+=1
        if r==rows: break
    return r
okg=True
for N in range(2,9):
    for s in range(0,N+1):
        for sp in range(s,N+1-s):
            A=list(itertools.combinations(range(N),s)); B=list(itertools.combinations(range(N),sp))
            M=[[1 if not set(a)&set(b) else 0 for b in B] for a in A]
            r=rank_exact(M) if len(A)*len(B)<=4000 else int(np.linalg.matrix_rank(np.array(M,dtype=float)))
            if r!=math.comb(N,min(s,sp)): okg=False; print("GOTTLIEB FAIL",N,s,sp,r)
print("Gottlieb: rank DISJ(s,s') = C(N,min(s,s')) for all N<=8, s+s'<=N:",okg)
# counter-check: s+s'>N gives zero matrix (so the hypothesis is needed)

# (c) Prop S
oks=True
for (n,d,K) in [(600,8,8),(600,8,30),(600,8,100),(900,16,16),(900,16,60),(900,16,300),(800,6,3),(800,6,40)]:
    G=nx.random_regular_graph(d,n,seed=rng.randrange(10**6))
    A=nx.to_numpy_array(G); ev=np.sort(np.abs(np.linalg.eigvalsh(A))); lam=ev[-2]
    for mode in ["random","bfs","unbalanced"]:
        if mode=="random":
            perm=list(range(n)); rng.shuffle(perm); c={v:i%K for i,v in enumerate(perm)}
        elif mode=="bfs":   # locality-preserving clustering: BFS order cut into K blocks (tries to minimise inc)
            order=list(nx.bfs_tree(G,0).nodes()); order+= [v for v in range(n) if v not in set(order)]
            sz=math.ceil(n/K); c={v:i//sz for i,v in enumerate(order)}
        else:               # K parts but one huge part
            perm=list(range(n)); rng.shuffle(perm); c={v:(i if i<K-1 else K-1) for i,v in enumerate(perm)}
        sizes=np.bincount([c[v] for v in range(n)]); smax=sizes.max()
        inc=sum(len({c[w] for w in G[u]}-{c[u]}) for u in range(n))
        # orientation: lower index -> higher index
        incp=sum(len({c[w] for w in G[u] if w>u}-{c[u]}) for u in range(n))
        lb=n*d*d/(d*d*smax/n+lam*lam)-n; lbp=n*d*d/(4*(d*d*smax/n+lam*lam))-n
        if inc<lb-1e-9 or incp<lbp-1e-9: oks=False
        print("PropS n=%d d=%d K=%d %-10s smax=%d lam=%.2f inc=%d lb=%.0f inc+=%d lb+=%.0f inc/(nd)=%.3f"%(n,d,K,mode,smax,lam,inc,lb,incp,lbp,inc/(n*d)))
print("Prop S bounds hold:",oks)
