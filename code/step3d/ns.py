import itertools, random, subprocess, sys, json, time
from math import comb
HERE='/tmp/claude-0/-home-user/2d12e629-e6d2-542b-bfc0-8b766f57383f/scratchpad/step3d'

def rand3cnf(n,m,seed):
    rng=random.Random(seed)
    return [tuple(v if rng.random()<.5 else -v for v in rng.sample(range(1,n+1),3)) for _ in range(m)]

def unsat(n,cl):
    # brute force over all 2^n assignments (bitmask)
    masks=[]
    for c in cl:
        # assignment falsifies c iff each literal false
        req=0; val=0
        for l in c:
            b=1<<(abs(l)-1); req|=b
            if l<0: val|=b   # x=1 falsifies -x
        masks.append((req,val))
    import numpy as np
    a=np.arange(1<<n,dtype=np.int64); alive=np.ones(1<<n,dtype=bool)
    for r,v in masks: alive &= (a&r)!=v
    return not alive.any()

def clause_poly(c):
    # p_C = prod over literals: x -> (1-x), -x -> x ; returns dict mask->coef
    p={0:1}
    for l in c:
        b=1<<(abs(l)-1); q={}
        for m,a in p.items():
            if l>0:
                q[m]=q.get(m,0)+a; q[m|b]=q.get(m|b,0)-a
            else:
                q[m|b]=q.get(m|b,0)+a
        p=q
    return {m:a for m,a in p.items() if a}

def ns_feasible(n,polys,d,q):
    # monomials of degree<=d, columns ordered by degree descending, constant last
    mons=[m for k in range(d,-1,-1) for S in itertools.combinations(range(n),k) for m in [sum(1<<i for i in S)]]
    idx={m:i for i,m in enumerate(mons)}
    rows=[]; trip=[]
    for p in polys:
        vs=0; dp=0
        for m in p: vs|=m
        dp=max(bin(m).count('1') for m in p)
        free=[i for i in range(n) if not vs>>i&1]
        for k in range(0,d-dp+1):
            for S in itertools.combinations(free,k):
                g=sum(1<<i for i in S); r=len(rows); rows.append(1)
                for m,a in p.items(): trip.append(f"{r} {idx[m|g]} {a%q}")
    inp=f"{q} {len(rows)} {len(mons)} {len(mons)-1} {len(trip)}\n"+"\n".join(trip)+"\n"
    out=subprocess.run([HERE+'/howell'],input=inp,capture_output=True,text=True).stdout
    return out.strip()=='1', len(rows), len(mons)

def min_deg(n,polys,q,dmax):
    dp=max(max(bin(m).count('1') for m in p) for p in polys)
    for d in range(dp,dmax+1):
        ok,_,_=ns_feasible(n,polys,d,q)
        if ok: return d
    return None  # > dmax

def tseitin(edges,charge):
    # vertex constraint XOR of incident edges = charge(v), as CNF
    V=sorted({u for e in edges for u in e}); cl=[]
    for v in V:
        inc=[i+1 for i,e in enumerate(edges) if v in e]
        for signs in itertools.product([0,1],repeat=len(inc)):
            # forbid assignment with parity != charge: assignment a_i=signs
            if sum(signs)%2!=charge.get(v,0):
                cl.append(tuple(-x if s else x for x,s in zip(inc,signs)))
    return cl
