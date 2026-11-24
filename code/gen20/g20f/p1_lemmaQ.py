# Lemma Q at scale: construct s (global max of Q), fractional x, LP vertex, round, check per-direction b_i >= M_i'/u - u
import numpy as np, itertools, random
from scipy.optimize import linprog
rng=np.random.default_rng(7); random.seed(7)
def run(n,M,u,r,skew):
    wires=[(random.randrange(n),random.randrange(n)) for _ in range(M)]
    labs=[np.array([min(r-1,int(rng.exponential(skew))) if skew else rng.integers(r) for _ in range(n)]) for _ in range(u)]
    # general r: sigma signs on groups with cut >= M_i/2 (best of all 2^r sign vectors)
    delta=np.zeros((M,u),dtype=np.int64); Mi=[]; Mip=[]
    for i,lab in enumerate(labs):
        cut=[(lab[h],lab[t]) for h,t in wires if lab[h]!=lab[t]]; Mi.append(len(cut))
        best=None
        for sg in itertools.product((-1,1),repeat=r):
            c=sum(sg[a]!=sg[b] for a,b in cut)
            if best is None or c>best[0]: best=(c,sg)
        sg=best[1]; Mip.append(best[0]); assert 2*best[0]>=len(cut)
        for e,(h,t) in enumerate(wires): delta[e,i]=(sg[lab[h]]-sg[lab[t]])//2
    S=np.array(list(itertools.product((-1,1),repeat=u)))
    Q=((delta@S.T)**2).sum(0); s=S[Q.argmax()]
    x=(delta@s)/u
    lhs=s*(x@delta)
    frac_ok=np.all(lhs>=np.array(Mip)/u-1e-9)
    # LP vertex of {x in [-1,1]^E : s_i sum_e x_e delta_e(i) >= Mip_i/u}  (simplex => vertex)
    res=linprog(c=rng.normal(size=M),A_ub=-(delta*s).T,b_ub=-np.array(Mip)/u,bounds=[(-1,1)]*M,method='highs-ds')
    xv=res.x; nfrac=np.sum((np.abs(np.abs(xv)-1)>1e-9))
    xi=np.where(np.abs(np.abs(xv)-1)<=1e-9,np.round(xv),0)
    # actual b_i: (1/2) sum_g |y(g)| with integer selection xi
    D=np.zeros((M,n)); 
    for e,(h,t) in enumerate(wires): D[e,h]+=1; D[e,t]-=1
    y=xi@D; b=[]
    for lab in labs:
        b.append(0.5*sum(abs(y[lab==g].sum()) for g in range(r)))
    b=np.array(b)
    return frac_ok, nfrac, np.all(b>=np.array(Mip)/u-u-1e-9), np.all(b>=np.array(Mi)/(2*u)-u-1e-9), (r==2 and not np.all(b>=np.array(Mi)/u-u-1e-9)), min(b-np.array(Mip)/u)
stats=dict(frac=0,nfrac_bad=0,rnd=0,gen=0,r2=0,N=0); worst=1e9
for trial in range(400):
    u=random.randint(1,8); r=random.choice([2,2,3,4,6]); n=random.randint(4,40); M=random.randint(5,200)
    f,nf,ok,okg,badr2,w=run(n,M,u,r,random.choice([0,0.5,1.5]))
    stats['N']+=1; stats['frac']+=not f; stats['nfrac_bad']+= nf>u; stats['rnd']+=not ok; stats['gen']+=not okg; stats['r2']+=badr2; worst=min(worst,w)
print(stats,' min_i (b_i - M_i\'/u) over all =',worst)
