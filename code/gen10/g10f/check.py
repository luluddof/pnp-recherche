import math, random
from collections import defaultdict
H=lambda x: -x*math.log(x)-(1-x)*math.log(1-x)
def kappa(beta,eps): return 3*(1-eps)*H(beta)-math.log(2)-H(eps)
for beta,eps in [(1/6,1/20),(1/8,1/50)]:
    k=kappa(beta,eps); print("beta",beta,"eps",eps,"H(beta)",H(beta),"H(eps)",H(eps),"kappa",k,"Delta>",math.log(2)/k)
def solve(rhs):
    lo,hi=1e-12,1.0
    for _ in range(200):
        mid=(lo+hi)/2
        if mid*math.log(math.e/mid)<rhs: lo=mid
        else: hi=mid
    return lo
for den in [9,60,63,150]:
    print("y log(e/y)=1/(%d e^2): y="%den, solve(1/(den*math.e**2)))
# deterministic sanity test of the transfer lemma
def test(n,m,seed):
    rnd=random.Random(seed)
    clauses=[tuple(rnd.sample(range(n),3)) for _ in range(m)]
    d=[0]*n
    copy={}
    Vs=[]; pi=[]
    adj=defaultdict(set)
    fib=defaultdict(list)
    for ei,e in enumerate(clauses):
        cs=[]
        for v in e:
            idx=len(Vs); Vs.append((v,d[v])); pi.append(v); fib[v].append(idx); d[v]+=1; cs.append(idx)
        for a in cs:
            for b in cs:
                if a!=b: adj[a].add(b)
    for v in range(n):
        f=fib[v]; rnd.shuffle(f)
        if len(f)>=2:
            for i in range(len(f)):
                a,b=f[i],f[(i+1)%len(f)]
                if a!=b: adj[a].add(b); adj[b].add(a)
    n0=0
    for v in range(n):
        if d[v]==0 and rnd.random()<0.5:
            idx=len(Vs); Vs.append((v,0)); pi.append(v); fib[v].append(idx); n0+=1
    ns=len(Vs)
    assert ns==3*m+n0
    # random Y*: BFS ball boundary
    for trial in range(20):
        k=rnd.randint(0,ns//6)
        Ys=set(rnd.sample(range(ns),k))
        # also add a BFS layer sometimes
        seen=set(Ys); comps=[]
        for s in range(ns):
            if s in seen: continue
            comp=[s]; seen.add(s); st=[s]
            while st:
                x=st.pop()
                for y in adj[x]:
                    if y not in seen: seen.add(y); comp.append(y); st.append(y)
            comps.append(comp)
        rnd.shuffle(comps)
        t=rnd.randint(0,len(comps))
        As=set(x for c in comps[:t] for x in c)
        Bs=set(range(ns))-As-Ys
        Y=set(pi[x] for x in Ys)
        A=set(); B=set()
        for v in range(n):
            if v in Y: continue
            f=fib[v]
            if f and all(x in As for x in f): A.add(v)
            else:
                assert all(x in Bs for x in f), "fibre split"
                B.add(v)
        L=sum(d[v] for v in Y)
        mA=mB=cross=0
        for e in clauses:
            if any(v in Y for v in e): 
                if any(v in A for v in e) and any(v not in A for v in e): cross+=1
                continue
            inA=all(v in A for v in e); inB=all(v in B for v in e)
            assert inA or inB
            mA+=inA; mB+=inB
        assert cross<=L
        assert 3*mA>=len(As)-3*L-n0-1e-9, (mA,len(As),L,n0)
        assert 3*mB>=len(Bs)-3*L-n0-1e-9
    return True
for s in range(200):
    test(random.Random(s).randint(8,40), random.Random(s+1000).randint(5,120), s)
print("transfer lemma sanity tests passed")
