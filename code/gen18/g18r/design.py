# Independent implementation of the per-point piecewise block design over GF(3^k).
import numpy as np, itertools, math
from gf import GF

def solve_hash(block, m):
    """a in F_3^m with sum_i a_i C(w,i) = (-1)^w mod 3 for all w in block. Brute force."""
    sols=[]
    for a in itertools.product(range(3), repeat=m):
        if all(sum(a[i]*math.comb(w,i) for i in range(m))%3 == ((-1)**w)%3 for w in block):
            sols.append(a)
    return sols

def blocks_for(Q, m, shift):
    """partition {0..Q} into intervals of length <= m; first block has length shift (1..m)."""
    bl=[]; lo=0; ln=shift
    while lo<=Q:
        hi=min(Q, lo+ln-1); bl.append(list(range(lo,hi+1))); lo=hi+1; ln=m
    return bl

class Design:
    def __init__(self, F, Q, m, shift, taus):
        self.F=F; self.Q=Q; self.m=m; self.taus=taus
        assert len(set(taus))==Q and 0 not in taus
        self.blocks=blocks_for(Q,m,shift)
        self.blk_of={w:i for i,b in enumerate(self.blocks) for w in b}
        self.coef=[]
        for b in self.blocks:
            s=solve_hash(b,m); assert s, ("(#) unsolvable", b)
            self.coef.append(s[0])
        self.tops=set(b[-1] for b in self.blocks if b[-1]!=Q)
        self._cache={}
    def Pz(self, z):
        """z: tuple of 0/1 (1 means sign -1). Returns polynomial in T (list over F), = (-1)^{|z|} F_beta(y(T))."""
        if z in self._cache: return self._cache[z]
        F=self.F; m=self.m
        # e_i(y), y_s = 1 - T/tau_s for s in supp z: gen function prod (1 + X y_s) truncated in X at m-1
        e=[[1]]+[[0] for _ in range(m-1)]   # e[i] = polynomial in T
        for s in range(self.Q):
            if z[s]:
                ys=[1, F.neg(F.inv(self.taus[s]))]
                for i in range(m-1,0,-1):
                    e[i]=F.padd(e[i], F.pmul(e[i-1], ys))
        a=self.coef[self.blk_of[sum(z)]]
        p=[0]
        for i in range(m):
            p=F.padd(p, F.pscal(F.emb(a[i]), e[i]))
        if sum(z)%2: p=F.pscal(F.neg(1), p)
        self._cache[z]=p
        return p
    def scalar_checks(self):
        """P_z(0)=1; deg<=m-1; pairing P_z(tau_t)+P_{z^t}(tau_t)=0 for good z_{-t}. Returns (maxdeg, #bad accidental good)"""
        F=self.F; Q=self.Q; maxdeg=-1; acc=0; nbadpairs=0
        for z in itertools.product((0,1), repeat=Q):
            p=self.Pz(z)
            assert F.peval(p,0)==1, ("P_z(0)!=1", z)
            maxdeg=max(maxdeg, F.pdeg(p))
            for t in range(Q):
                if z[t]==0:
                    z2=z[:t]+(1,)+z[t+1:]
                    s=F.add(F.peval(p,self.taus[t]), F.peval(self.Pz(z2), self.taus[t]))
                    w=sum(z)
                    if w not in self.tops:
                        assert s==0, ("pairing fails", z, t)
                    else:
                        nbadpairs+=1
                        if s==0: acc+=1
        assert maxdeg<=self.m-1
        return maxdeg, nbadpairs, acc
    def bad_fraction(self, r):
        """Pr[W' in tops], W' ~ Bin(Q-1,(r-1)/2r), exact rational"""
        from fractions import Fraction
        p=Fraction(r-1,2*r)
        return sum(math.comb(self.Q-1,w)*p**w*(1-p)**(self.Q-1-w) for w in self.tops if w<=self.Q-1)

def eigvecs(r):
    """list of (vector on [r]^2 as int array length r*r (index a*r+d), sign 0/1, normsq)"""
    out=[]
    for a in range(r):
        v=np.zeros(r*r,dtype=np.int64); v[a*r+a]=1; out.append((v,0,1))
    for a in range(r):
        for d in range(a+1,r):
            v=np.zeros(r*r,dtype=np.int64); v[a*r+d]=1; v[d*r+a]=1; out.append((v,0,2))
            v=np.zeros(r*r,dtype=np.int64); v[a*r+d]=1; v[d*r+a]=-1; out.append((v,1,2))
    return out
