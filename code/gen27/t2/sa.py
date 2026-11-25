# simulated annealing in the MC model (signed partial permutations): state g_l, h_l (maps [N] -> [N] u {-1}), signs.
# Theta[(x,y),(h(y),g(x))] += sg(x) sh(y); cost = rank_{GF(P)}(1 - Theta).  usage: sa.py N q steps seed [T0] [init.json]
import os
os.environ['OMP_NUM_THREADS']='1'
import sys, numpy as np, flint, math, time, json
P=1000003
def theta(N,g,h,sg,sh):
    n=N*N; T=np.zeros((n,n),dtype=np.int64)
    for l in range(len(g)):
        xs=np.where(g[l]>=0)[0]; ys=np.where(h[l]>=0)[0]
        if len(xs)==0 or len(ys)==0: continue
        X,Y=np.meshgrid(xs,ys,indexing='ij')
        np.add.at(T,(X*N+Y,h[l][Y]*N+g[l][X]),sg[l][X]*sh[l][Y])
    return T
def cost(N,st):
    T=theta(N,*st); M=np.eye(N*N,dtype=np.int64)-T
    return flint.nmod_mat((M%P).tolist(),P).rank()
def rand_pp(N,rng,dens):
    perm=rng.permutation(N); m=np.where(rng.random(N)<dens,perm,-1); return m
def mutate(N,st,rng):
    g,h,sg,sh=[ [a.copy() for a in s] for s in st]
    l=rng.integers(len(g)); side=rng.integers(2); m=(g,h)[side][l]; sgn=(sg,sh)[side][l]
    r=rng.random(); x=rng.integers(N)
    if r<0.2: sgn[x]*=-1
    elif r<0.8:
        y=rng.integers(N); xs=np.where(m==y)[0]
        if len(xs): m[xs[0]]=m[x]
        m[x]=y
    elif r<0.9: m[x]=-1
    else:
        free=[y for y in range(N) if y not in m]
        if free and m[x]<0: m[x]=rng.choice(free)
        else: sgn[x]*=-1
    return (g,h,sg,sh)
def run(N,q,steps,seed,T0=0.6,init=None,log=True):
    rng=np.random.default_rng(seed)
    if init is not None: st=init
    else: st=([rand_pp(N,rng,1.0) for _ in range(q)],[rand_pp(N,rng,1.0) for _ in range(q)],[rng.choice([-1,1],N) for _ in range(q)],[rng.choice([-1,1],N) for _ in range(q)])
    c=cost(N,st); best=(c,st); t0=time.time()
    for it in range(steps):
        T=T0*(1-it/steps)+0.02
        nst=mutate(N,st,rng); nc=cost(N,nst)
        if nc<=c or rng.random()<math.exp(-(nc-c)/T):
            st,c=nst,nc
            if c<best[0]:
                best=(c,st)
                if log: print(f"  it {it} rank {c} t={time.time()-t0:.0f}s",flush=True)
    return best
def to_EC(N,st):
    g,h,sg,sh=st; E=[];C=[]
    for l in range(len(g)):
        e=np.zeros((N,N),dtype=int); c=np.zeros((N,N),dtype=int)
        for x in range(N):
            if g[l][x]>=0: e[x,g[l][x]]=sg[l][x]
            if h[l][x]>=0: c[x,h[l][x]]=sh[l][x]
        E.append(e); C.append(c)
    return E,C
if __name__=='__main__':
    N,q,steps,seed=map(int,sys.argv[1:5]); T0=float(sys.argv[5]) if len(sys.argv)>5 else 0.6
    b=run(N,q,steps,seed,T0)
    print("BEST",b[0])
    from lib import save, pq
    E,C=to_EC(N,b[1]); p,qq=pq(N,E,C,exact=True); print("exact check",p,qq); save(f"work/sa_{N}_{p}_{qq}_s{seed}.json",N,E,C)
