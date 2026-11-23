# Threshold-of-local-oracles distinguishers: y -> sign(sum_S w_S[pattern of y on S]), |S|<=q
# This is exactly the "circuit with small fan-in oracle gates + one top threshold gate" shape from the locality barrier.
import json, itertools, numpy as np, scipy.sparse as sp, time
from sklearn.linear_model import LogisticRegression
n=4; N=16
c4=np.array(json.load(open('c4.json')))
F=np.arange(1<<N); B=((F[:,None]>>np.arange(N))&1)
rng=np.random.default_rng(0); perm=rng.permutation(1<<N)
tr,te=perm[:40000],perm[40000:]
def design(q, idx):
    subs=list(itertools.combinations(range(N),q))
    rows=[];cols=[]
    Bi=B[idx]
    for si,S in enumerate(subs):
        pat=(Bi[:,list(S)]*(1<<np.arange(q))).sum(1)
        cols.append(si*(1<<q)+pat)
    cols=np.stack(cols,1).ravel()
    rows=np.repeat(np.arange(len(idx)),len(subs))
    return sp.csr_matrix((np.ones(len(cols),dtype=np.float32),(rows,cols)),shape=(len(idx),len(subs)<<q))
for k in (4,5):
    y=(c4<=k)
    print(f"target C<={k} (base rate {y.mean():.3f})")
    for q in (1,2,3,4):
        t0=time.time()
        Xtr,Xte=design(q,tr),design(q,te)
        best=0
        for Cr in (0.1,1.0):
            m=LogisticRegression(C=Cr,max_iter=300,class_weight='balanced').fit(Xtr,y[tr])
            p=m.predict(Xte); yt=y[te]
            bal=0.5*((p&yt).sum()/yt.sum()+((~p)&(~yt)).sum()/(~yt).sum())
            best=max(best,bal)
        print(f"  q={q}: #oracle gates={Xtr.shape[1]>>q}, test balanced acc={best:.3f} ({time.time()-t0:.0f}s)")
