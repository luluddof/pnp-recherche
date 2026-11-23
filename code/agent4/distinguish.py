import json, itertools, numpy as np
from collections import Counter
n=4; N=16
c4=np.array(json.load(open('c4.json'))); c3=json.load(open('c3.json'))
F=np.arange(1<<N); B=((F[:,None]>>np.arange(N))&1)
# ---------- restriction (subcube) test: C(f|rho) <= C(f) always ----------
def restrict_tt(f, var, val):
    out=0; j=0
    for t in range(N):
        if ((t>>var)&1)==val:
            out|=((f>>t)&1)<<j; j+=1
    return out
R1=np.zeros(1<<N,dtype=int)
for f in range(1<<N):
    R1[f]=max(c3[restrict_tt(f,v,b)] for v in range(n) for b in (0,1))
print("== Restriction test (max complexity over 8 one-variable restrictions, uses n=3 MCSP oracle on 8-bit windows) ==")
for k in range(2,7):
    yes=c4<=k; no=~yes
    passed=R1<=k
    fp=(passed&no).sum()/no.sum()
    print(f"k={k}: P(yes)={yes.mean():.3f}  NO-instances wrongly passing local test: {fp:.3f}")
print("joint (C4, maxC3restr):",sorted(Counter(zip(c4.tolist(),R1.tolist())).items()))
# lower bound C4 >= maxR + 1 when f depends on restricted var? check C4 - R1 distribution
print("C4 - maxR1 distribution:",sorted(Counter((c4-R1).tolist()).items()))
# ---------- global features ----------
pm=1-2*B
H=np.array([[(-1)**bin(s&t).count('1') for t in range(N)] for s in range(N)])
four=pm@H.T/N
feats={}
feats['weight']=B.sum(1)
feats['balance |w-8|']=np.abs(B.sum(1)-8)
feats['total influence']=sum((B!=B[:,[t^(1<<i) for t in range(N)]]).sum(1) for i in range(n))/N
feats['spectral L1']=np.abs(four).sum(1)
# GF(2) degree via Moebius
A=B.copy()
for i in range(n):
    for t in range(N):
        if (t>>i)&1: A[:,t]^=A[:,t^(1<<i)]
deg=np.zeros(1<<N,dtype=int)
for t in range(N):
    deg=np.maximum(deg,A[:,t]*bin(t).count('1'))
feats['GF2 degree']=deg
feats['ANF #monomials']=A.sum(1)
# distinct subfunctions (Nechiporuk/communication-matrix rank-ish), max over 6 row/col splits
ds=np.zeros(1<<N,dtype=int)
for rv in itertools.combinations(range(n),2):
    cv=[i for i in range(n) if i not in rv]
    rows=[]
    for a in range(4):
        idx=[]
        for b in range(4):
            t=0
            t|=(a&1)<<rv[0]; t|=((a>>1)&1)<<rv[1]; t|=(b&1)<<cv[0]; t|=((b>>1)&1)<<cv[1]
            idx.append(t)
        rows.append(B[:,idx]@np.array([1,2,4,8]))
    rows=np.stack(rows,1)
    nd=np.array([len(set(r)) for r in rows.tolist()])
    ds=np.maximum(ds,nd)
feats['max #distinct 2-var subfns']=ds
feats['maxC over 1-var restrictions']=R1
def best_bal_acc(x,y):
    best=0
    for th in np.unique(x):
        for pred in (x<=th, x>th):
            acc=0.5*((pred&y).sum()/y.sum()+((~pred)&(~y)).sum()/(~y).sum())
            best=max(best,acc)
    return best
from scipy.stats import spearmanr
y=c4<=5
print("\n== Single-feature distinguishers, target C<=5 vs C>=6 (balanced accuracy; 0.5=chance) ==")
for k,v in feats.items():
    print(f"{k:32s} bal.acc={best_bal_acc(v,y):.3f}  spearman(C)={spearmanr(v,c4)[0]:+.3f}")
np.save('feats.npy',np.stack([feats[k] for k in feats],1))
