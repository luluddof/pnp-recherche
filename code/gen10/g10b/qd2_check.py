# Independent exhaustive check of the G8d Q_{d+2} lemma (own code).
# Unknowns c_{z,i} (i=1..d), P_z = 1 + sum_i c_{z,i} t^i.  Edge (z,z^k) good <=> 2 + sum_i (c_{z,i}+c_{z^k,i}) tau_k^i = 0.
# Good set G feasible <=> affine system consistent.  Lemma: no G = (all edges) minus 3 edges is feasible.
import itertools, sys, random
from fractions import Fraction as Fr
P=(1<<61)-1

def consistent(rows, mod):
    # rows: list of (coeff list, rhs). Gaussian elimination; returns True iff consistent.
    M=[r[:] for r in rows]; nc=len(M[0])-1; r=0
    for c in range(nc):
        p=None
        for i in range(r,len(M)):
            if M[i][c]!=0: p=i;break
        if p is None: continue
        M[r],M[p]=M[p],M[r]
        if mod: inv=pow(M[r][c],mod-2,mod)
        else: inv=1/M[r][c]
        for i in range(r+1,len(M)):
            if M[i][c]!=0:
                f=M[i][c]*inv
                if mod: M[i]=[(a-f*b)%mod for a,b in zip(M[i],M[r])]
                else: M[i]=[a-f*b for a,b in zip(M[i],M[r])]
        r+=1
        if r==len(M): break
    for i in range(r,len(M)):
        if M[i][nc]!=0: return False
    return True

def run(d,tau,mod,ks=(3,)):
    l=d+2; n=1<<l
    edges=[(z,k) for k in range(l) for z in range(n) if not (z>>k)&1]
    def row(z,k):
        r=[0]*(n*d+1)
        for i in range(1,d+1):
            v=pow(tau[k],i,mod) if mod else tau[k]**i
            r[z*d+i-1]=v; r[(z^(1<<k))*d+i-1]=v
        r[n*d]=(-2)%mod if mod else Fr(-2)
        return r
    rows={e:row(*e) for e in edges}
    idx={e:i for i,e in enumerate(edges)}
    res={}
    for nb in ks:
        feas=0; tot=0; seen=set()
        for B in itertools.combinations(edges,nb):
            # translation symmetry z -> z xor a : canonical rep
            key=min(tuple(sorted(idx[((z^a)&~(1<<k)),k] for z,k in B)) for a in range(n))
            if key in seen: continue
            seen.add(key); tot+=1
            Bs=set(B)
            if consistent([rows[e] for e in edges if e not in Bs],mod): feas+=1
        res[nb]=(feas,tot)
    return res

if __name__=="__main__":
    random.seed(5)
    # d=0: all 4 edges of Q_2 bad trivially
    for d,taus in [(1,[[Fr(1),Fr(2),Fr(3)],[Fr(1),Fr(-1),Fr(2)],[Fr(1),Fr(1,2),Fr(1,3)],[Fr(3,7),Fr(-5,2),Fr(11)]]),
                   (2,[[Fr(1),Fr(2),Fr(3),Fr(4)],[Fr(1),Fr(-1),Fr(2),Fr(-2)],[Fr(1),Fr(1,2),Fr(1,3),Fr(1,4)]])]:
        for tau in taus:
            r=run(d,tau,0,ks=(0,1,2,3,4))
            print("d=%d tau=%s exact/Q  {#bad: (feasible orbits, orbits)} = %s"%(d,[str(t) for t in tau],r)); sys.stdout.flush()
    for tau in [[1,2,3,4,5],[1,P-1,2,P-2,3],[random.randrange(1,P) for _ in range(5)]]:
        r=run(3,tau,P,ks=(1,2,3))
        print("d=3 tau=%s mod p {#bad: (feasible orbits, orbits)} = %s"%(tau,r)); sys.stdout.flush()
