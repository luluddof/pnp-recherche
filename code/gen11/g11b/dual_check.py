"""G11b: checks of the dual (Smolensky-type) lower bound.
Commutative class: P_z(t), deg<=d, P_z(0)=1; k-edge {z,z^k} good iff (P_z+P_{z^k})(tau_k)=0.
Certificate space W_m = { w(z,t) = sum_{|T|<=m} s_T prod_{j in T}(1-t/tau_j) r_T(t), deg r_T <= m-|T| }.
Claim: feasible  =>  evaluation map W_m -> R^{bad edges}, w -> (w(z,tau_k)), is injective, m=floor((l-d-1)/2).
Test 1: layered (G7a) feasible configurations: rank of evaluation map == dim W_m (mod p).
Test 2: bad set with nontrivial kernel -> explicit Fredholm certificate y (exact over Q) and primal infeasible.
Test 3: general matrix class: identity a(tau_k) K_k (aX)(tau_k)^T = 0 for weakly nice K_k (exact over Q)."""
import itertools, random, sys
from fractions import Fraction as Fr
from math import comb
P=(1<<31)-1
def rank_mod(rows, p=P):
    M=[r[:] for r in rows]; r=0; nc=len(M[0]) if M else 0
    for c in range(nc):
        piv=None
        for i in range(r,len(M)):
            if M[i][c]%p: piv=i;break
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]; inv=pow(M[r][c],p-2,p); row=[x*inv%p for x in M[r]]; M[r]=row
        for i in range(len(M)):
            if i!=r and M[i][c]%p:
                f=M[i][c]; M[i]=[(a-f*b)%p for a,b in zip(M[i],row)]
        r+=1
        if r==len(M): break
    return r
def basis(l,m):
    return [(T,i) for s in range(m+1) for T in itertools.combinations(range(l),s) for i in range(m-s+1)]
def evalrow(z,k,l,tau,B,field):
    # value at (z,tau_k) of basis element s_T prod_{j in T}(1-t/tau_j) t^i
    row=[]
    for T,i in B:
        if k in T: row.append(0); continue
        v=tau[k]**i
        for j in T:
            sj=-1 if (z>>j)&1 else 1
            if field=='Q': v*=sj*(1-Fr(tau[k])/tau[j])
            else: v=v*sj*(1-tau[k]*pow(tau[j],P-2,P))%P
        row.append(v if field=='Q' else v%P)
    return row
def layered_bad(l,d,r):
    cuts=set(); c=r
    while c<l: cuts.add(c); c+=d+1      # boundary between weight c and c+1
    return [(z,k) for z in range(1<<l) for k in range(l) if not (z>>k)&1 and bin(z).count('1') in cuts]
def Nm(l,d):
    m=(l-d-1)//2
    return sum(sum(comb(l,i) for i in range(j+1)) for j in range(m+1)) if m>=0 else 0
def test1(lmax=10):
    print("Test 1: layered configurations, rank of evaluation map vs dim W_m")
    for l in range(3,lmax+1):
        for d in range(1,min(l,6)):
            m=(l-d-1)//2
            if m<0: continue
            B=basis(l,m); dim=len(B); assert dim==Nm(l,d)
            best=min((len(layered_bad(l,d,r)),r) for r in range(d+1))
            bad=layered_bad(l,d,best[1])
            if dim>1200: print(f" l={l} d={d} m={m} dim={dim} #bad={len(bad)} (rank skipped)"); continue
            tau=[random.randrange(2,P) for _ in range(l)]
            # use a random subset of rows of size min(#bad, dim+40) -- rank must reach dim on ALL bad rows; subsample only for speed if huge
            rows=[evalrow(z,k,l,tau,B,'p') for z,k in bad]
            rk=rank_mod(rows)
            print(f" l={l} d={d} m={m} dim W_m={dim} #bad(layered)={len(bad)} rank={rk} injective={rk==dim}",flush=True)
            assert rk==dim
def solve_kernel_Q(rows,nc):
    M=[r[:] for r in rows]; piv=[]; r=0
    for c in range(nc):
        p=None
        for i in range(r,len(M)):
            if M[i][c]!=0: p=i;break
        if p is None: continue
        M[r],M[p]=M[p],M[r]; inv=1/M[r][c]; M[r]=[x*inv for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]!=0:
                f=M[i][c]; M[i]=[a-f*b for a,b in zip(M[i],M[r])]
        piv.append(c); r+=1
        if r==len(M): break
    free=[c for c in range(nc) if c not in piv]
    if not free: return None
    f=free[0]; v=[Fr(0)]*nc; v[f]=Fr(1)
    for i,c in enumerate(piv): v[c]=-M[i][f]
    return v
def primal_feasible(l,d,tau,badset):
    n=1<<l; rows=[]
    for z in range(n):
        for k in range(l):
            if (z>>k)&1 or (z,k) in badset: continue
            row=[0]*(n*d+1)
            for i in range(1,d+1):
                v=pow(tau[k],i,P); row[z*d+i-1]=v; row[(z|1<<k)*d+i-1]=v
            row[n*d]=P-2; rows.append(row)
    A=[r[:-1] for r in rows]
    return rank_mod(A)==rank_mod(rows)
def test2(l,d,trials=3,seed=1):
    random.seed(seed); m=(l-d-1)//2; B=basis(l,m); dim=len(B); tau=[Fr(k+1) for k in range(l)]
    edges=[(z,k) for z in range(1<<l) for k in range(l) if not (z>>k)&1]
    print(f"Test 2: (l,d)=({l},{d}) m={m} dim W_m={dim}: random bad sets of size dim-1 -> explicit certificate")
    for tr in range(trials):
        bad=random.sample(edges,dim-1); badset=set(bad)
        rows=[evalrow(z,k,l,tau,B,'Q') for z,k in bad]
        v=solve_kernel_Q(rows,dim); assert v is not None
        def w(z,t):
            tot=Fr(0)
            for c,(T,i) in zip(v,B):
                if c==0: continue
                x=c*t**i
                for j in T: x*=(-1 if (z>>j)&1 else 1)*(1-t/tau[j])
                tot+=x
            return tot
        # strip powers of t if w(.,0)==0 (did not happen generically); certificate W = w^2
        s0=sum(w(z,Fr(0))**2 for z in range(1<<l)); assert s0>0, "w(.,0)=0: divide by t first"
        mu=[]
        for k in range(l):
            x=tau[k]
            for j in range(l):
                if j!=k: x*=(tau[k]-tau[j])
            mu.append(1/x)
        ok=True; tot=Fr(0)
        for z in range(1<<l):
            ys=[]
            for k in range(l):
                a=w(z,tau[k]); b=w(z^(1<<k),tau[k])
                if a!=b: ok=False
                e=(z&~(1<<k),k)
                if e in badset and a!=0: ok=False
                ys.append(mu[k]*a*a)
            for i in range(1,d+1):
                if sum(y*tau[k]**i for k,y in enumerate(ys))!=0: ok=False
            tot+=sum(ys)
        tot/=2
        feas=primal_feasible(l,d,[int(t) for t in tau],badset)
        print(f"  trial {tr}: certificate valid (edge-consistent, zero on bad, vertex moments 1..d vanish)={ok}, sum y={float(tot):.3e} (nonzero={tot!=0}), primal feasible(mod p)={feas}",flush=True)
        assert ok and tot!=0 and not feas
def test3(l,d,trials=20,seed=2):
    random.seed(seed); n=1<<l; m=(l-d-1)//2; tau=[Fr(k+1) for k in range(l)]
    okall=True
    for tr in range(trials):
        Jc=[[[Fr(random.randint(-3,3)) for _ in range(n)] for _ in range(n)] for _ in range(d+1)]
        for x in range(n):
            for y in range(n): Jc[0][x][y]=Fr(1 if y==(n-1)^x else 0)
        r={x:[Fr(random.randint(-3,3)) for _ in range(m-bin(x).count('1')+1)] for x in range(n) if bin(x).count('1')<=m}
        def a(x,t):
            if x not in r: return Fr(0)
            v=sum(c*t**i for i,c in enumerate(r[x]))
            for k in range(l):
                if (x>>k)&1: v*=(1-t/tau[k])
            return v
        def J(t): return [[sum(Jc[i][x][y]*t**i for i in range(d+1)) for y in range(n)] for x in range(n)]
        def Phi(t,Mt): return sum(a(x,t)*Mt[x][y]*a((n-1)^y,t) for x in range(n) for y in range(n))
        if Phi(Fr(0),J(Fr(0)))!=sum(a(x,Fr(0))**2 for x in range(n)): okall=False
        for k in range(l):
            Jt=J(tau[k])
            K=[[Jt[x][y] if ((x^y)>>k)&1==0 else Fr(0) for y in range(n)] for x in range(n)]
            E=[[Jt[x][y]-K[x][y] for y in range(n)] for x in range(n)]
            if Phi(tau[k],K)!=0: okall=False
            if Phi(tau[k],Jt)!=Phi(tau[k],E): okall=False
    print(f"Test 3: general matrix class (l={l},d={d},m={m}): a(tau_k) K_k (aX)^T = 0 and Phi(0)=|a(0)|^2 in {trials} random trials: {okall}")
    assert okall
if __name__=="__main__":
    random.seed(0)
    test3(3,0); test3(4,1); test3(3,2)
    test2(6,1); test2(7,2,trials=2); test2(6,3,trials=3)
    test1(int(sys.argv[1]) if len(sys.argv)>1 else 9)
