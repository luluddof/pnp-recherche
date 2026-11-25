# Tabu search over uniform-type supports S (|S|=q): maximise #{w: f_S(w)=c}, f_S(w)=sum_{s in S} (-1)^{Q(s+w)}; exact ints.
# Optionally weights in {+1,-1} per element (signed design): f = sum sgn_s (-1)^{Q(s+w)}.
import numpy as np, sys, time
from bell import beta
k=int(sys.argv[1]); q=int(sys.argv[2]); steps=int(sys.argv[3]); seed=int(sys.argv[4]); signed=int(sys.argv[5]) if len(sys.argv)>5 else 0
rng=np.random.default_rng(seed); n=4**k; bt=beta(k).astype(np.int16); z=np.arange(n)
Bm=bt[z[:,None]^z[None,:]]          # Bm[s,w] = beta(s+w)
cand=np.concatenate([Bm,-Bm]) if signed else Bm   # candidate rows (signed: index n+s means -beta(s+.))
nc=len(cand)
S=list(rng.choice(nc,q,replace=False)); f=cand[S].sum(0).astype(np.int32)
def modecount(F):   # F: (m,n) ints in [-q,q]; best nonzero-value count per row
    m=F.shape[0]; off=(F+q)+(2*q+1)*np.arange(m)[:,None]
    bc=np.bincount(off.ravel(),minlength=m*(2*q+1)).reshape(m,2*q+1); bc[:,q]=-1
    return bc.max(1), bc.argmax(1)-q
cur=modecount(f[None])[0][0]; best=(cur,sorted(S)); tabu={}
t0=time.time()
for st in range(steps):
    bestmove=None
    order=rng.permutation(q)
    for i in order[:min(q,8)]:
        F=(f-cand[S[i]])[None,:]+cand
        sc,_=modecount(F); sc=sc.astype(np.int64)
        sc[S]=-10**9
        for t,v in tabu.items():
            if v>st and t<nc: sc[t]=-10**9
        j=int(np.argmax(sc + rng.random(nc)*0.5))
        if bestmove is None or sc[j]>bestmove[0]: bestmove=(sc[j],i,j)
    sc,i,j=bestmove
    tabu[S[i]]=st+max(5,q//2); f=f-cand[S[i]]+cand[j]; S[i]=j; cur=sc
    if cur>best[0]: best=(cur,sorted(S))
print(f'k={k} q={q} signed={signed} best agreement {best[0]}/{n} p={n-best[0]} ({time.time()-t0:.0f}s) S={[int(x) for x in best[1]]}',flush=True)
