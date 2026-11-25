# sanity checks: (1) Hajek identity of R3 and martingale drift; (2) stratified vs plain random grouping (R2);
# (3) time-space line of R4 on the discrete clique (outer set S = middle square, sweep rest).
import numpy as np, math, itertools
rng = np.random.default_rng(1)
def net(g, kind):
    if kind=='unif': c = rng.random((g,g))
    elif kind=='hetero':   # heterogeneous node sums: weight ~ u_a u_b
        u = rng.random(g)**3; c = np.outer(u,u)*rng.random((g,g))
    elif kind=='block':
        lab = rng.integers(0,2,g); c = np.where(lab[:,None]==lab[None,:], 1.0, 0.1)*rng.random((g,g))
    c = np.triu(c,1); c = c+c.T; return c/c.max()      # pair weights <= 1 = P/g
# (1)
g=40; c=net(g,'hetero'); U=c.sum(); mu=U/(g*(g-1)); ca=c.sum(1)
al=(ca-(g-1)*mu)/(g-2); ga=c-mu-al[:,None]-al[None,:]; np.fill_diagonal(ga,0)
print("row sums of gamma ~0:", np.abs(ga.sum(1)).max())
perm=rng.permutation(g); err=0
for j in range(1,g):
    S=perm[:j]; T=perm[j:]
    F=c[np.ix_(S,T)].sum(); A=al[S].sum(); Gm=ga[np.ix_(S,T)].sum()
    err=max(err,abs(F-(j*(g-j)*mu+(g-2*j)*A+Gm)))
print("Hajek identity max error:", err)
# drift E[Gamma_{j+1}|S_j] = Gamma_j (1-2/(g-j))
S=list(perm[:15]); T=list(perm[15:]); G0=ga[np.ix_(S,T)].sum()
nxt=[ga[np.ix_(S+[v],[w for w in T if w!=v])].sum() for v in T]
print("drift check:", np.mean(nxt), G0*(1-2/(g-15)))
# max cut of random orders vs U/4
for kind in ('unif','hetero','block'):
    g=512; c=net(g,kind); U=c.sum(); best=[]
    for r in range(5):
        p=rng.permutation(g); M=c[np.ix_(p,p)]
        # prefix cuts: F_j = F_{j-1} + c_v - 2 c(v, prefix)
        cs=M.sum(1); pre=np.cumsum(np.tril(M,-1).sum(1)); F=np.cumsum(cs)-2*pre
        best.append(F.max()/(U/4))
    print(f"order lemma g={g} {kind}: max cut/(U/4) over 5 random orders: {np.round(best,4)}  (P/U*4 = {4*c.sum(1).max()/U:.4f})")
# (2) grouping
for kind in ('unif','hetero'):
    K=16; t=64; g=K*t; c=net(g,kind); U=c.sum(); ca=c.sum(1)
    def maxY(groups):
        Gm=np.zeros((g,K)); Gm[np.arange(g),groups]=1
        Y=Gm.T@c@Gm; np.fill_diagonal(Y,-1); return Y.max()*K*K/U-1
    plain=[maxY(rng.permutation(np.repeat(np.arange(K),t))) for _ in range(20)]
    order=np.argsort(ca); strat=[]
    for _ in range(20):
        gr=np.empty(g,int)
        for i in range(t): gr[order[i*K:(i+1)*K]]=rng.permutation(K)
        strat.append(maxY(gr))
    print(f"grouping {kind} K={K} t={t}: (K^2 maxY/U - 1) plain mean {np.mean(plain):.4f}, stratified mean {np.mean(strat):.4f}")
# (3) R4: discrete time/memory of 'outer square + sweep' on K-clique
def tm(K, d):
    # nodes 0..K-1, S = edges {u,v}: m-d <= u < m <= v < m+d ; sweep in order 0..K-1
    m=K//2; S={(u,v) for u in range(m-d,m) for v in range(m,m+d)}
    T=0; Mx=0
    for j in range(K):   # step j: enumerate cut C_j (u<j<=v... ) plus new edges at node j
        cut={(u,v) for u in range(j) for v in range(j,K)}-S
        new={(j,v) for v in range(j+1,K)}-S
        T=max(T,len(cut|new)); Mx=max(Mx,len({(u,v) for u in range(j+1) for v in range(j+1,K)}-S))
    return len(S)+T, Mx
K=40
for d in (0,4,8,12,16,20):
    t,mm=tm(K,d); print(f"R4 K={K} d={d}: time {t/K**2:.4f} mem {mm/K**2:.4f}  line (1-M)/3 = {(1-mm/K**2)/3:.4f}")
