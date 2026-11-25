# Simulated annealing over sets S (|S|=q) maximising #{w : f_S(w) = c} (best c != 0), f_S = 1_S * beta.  Exact integers.
import numpy as np, sys, time
from bell import beta
k=int(sys.argv[1]); q=int(sys.argv[2]); iters=int(sys.argv[3]); seed=int(sys.argv[4]) if len(sys.argv)>4 else 0
T0=float(sys.argv[5]) if len(sys.argv)>5 else 2.0
rng=np.random.default_rng(seed); n=4**k; bt=beta(k).astype(np.int32)
z=np.arange(n)
def shift(s): return bt[z ^ s]           # beta(s + w) as function of w
S=list(rng.choice(n,q,replace=False)); inS=np.zeros(n,bool); inS[S]=True
f=np.zeros(n,np.int32)
for s in S: f+=shift(s)
def score(f):
    vals,cnt=np.unique(f,return_counts=True); m=vals!=0
    i=np.argmax(np.where(m,cnt,-1)); return int(cnt[i]), int(vals[i])
cur,cv=score(f); best=(cur,cv,sorted(S)); t0=time.time()
for it in range(iters):
    T=T0*(1-it/iters)+0.05
    i=rng.integers(q); s_old=S[i]; s_new=int(rng.integers(n))
    if inS[s_new]: continue
    f2=f-shift(s_old)+shift(s_new); sc,c2=score(f2)
    if sc>=cur or rng.random()<np.exp((sc-cur)/T*(64.0/n)*16):
        f=f2; S[i]=s_new; inS[s_old]=False; inS[s_new]=True; cur,cv=sc,c2
        if cur>best[0]: best=(cur,cv,sorted(S))
print(f'k={k} q={q} best agreement {best[0]}/{n}  p={n-best[0]}  c={best[1]}  ({time.time()-t0:.0f}s) S={best[2]}')
