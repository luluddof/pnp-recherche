# GF(3^k) via tables; elements are ints 0..3^k-1 (base-3 digits = coefficients)
import numpy as np, itertools
class GF:
    def __init__(self, k):
        # modulus: x^2+1 (k=2), x^3-x-1 i.e. x^3 = x+1 (k=3), k=1 trivial
        self.k=k; self.q=3**k
        red={1:[0],2:[2,0],3:[1,1,0]}[k]   # x^k = sum red[i] x^i
        def dig(a): return [(a//3**i)%3 for i in range(k)]
        def und(d): return sum((d[i]%3)*3**i for i in range(k))
        q=self.q
        self.ADD=np.zeros((q,q),dtype=np.int64); self.MUL=np.zeros((q,q),dtype=np.int64)
        for a in range(q):
            for b in range(q):
                da,db=dig(a),dig(b)
                self.ADD[a,b]=und([da[i]+db[i] for i in range(k)])
                prod=[0]*(2*k)
                for i in range(k):
                    for j in range(k): prod[i+j]+=da[i]*db[j]
                for e in range(2*k-1,k-1,-1):
                    c=prod[e]%3; prod[e]=0
                    for i in range(k): prod[e-k+i]+=c*(red[i] if i<len(red) else 0)
                self.MUL[a,b]=und(prod[:k])
        self.NEG=np.array([und([-x for x in dig(a)]) for a in range(q)])
        self.INV=np.zeros(q,dtype=np.int64)
        for a in range(1,q):
            for b in range(1,q):
                if self.MUL[a,b]==1: self.INV[a]=b
        assert all(self.INV[a]>0 for a in range(1,q))
        # checks: field axioms sample
        assert all(self.MUL[a,self.INV[a]]==1 for a in range(1,q))
    def add(s,a,b): return int(s.ADD[a,b])
    def mul(s,a,b): return int(s.MUL[a,b])
    def neg(s,a): return int(s.NEG[a])
    def sub(s,a,b): return int(s.ADD[a,s.NEG[b]])
    def inv(s,a): assert a!=0; return int(s.INV[a])
    def emb(s,n): return n%3          # integer -> prime field
    # polynomials: lists of elements, low degree first
    def padd(s,p,q_):
        n=max(len(p),len(q_)); p=p+[0]*(n-len(p)); q_=q_+[0]*(n-len(q_))
        return [s.add(a,b) for a,b in zip(p,q_)]
    def pmul(s,p,q_,trunc=None):
        n=len(p)+len(q_)-1
        if trunc: n=min(n,trunc)
        out=[0]*n
        for i,a in enumerate(p):
            if a==0: continue
            for j,b in enumerate(q_):
                if i+j>=n: break
                out[i+j]=s.add(out[i+j],s.mul(a,b))
        return out
    def pscal(s,c,p): return [s.mul(c,a) for a in p]
    def peval(s,p,x):
        r=0
        for a in reversed(p): r=s.add(s.mul(r,x),a)
        return r
    def pdeg(s,p):
        d=-1
        for i,a in enumerate(p):
            if a: d=i
        return d
if __name__=="__main__":
    for k in (1,2,3):
        F=GF(k); q=F.q
        # associativity/distributivity exhaustive
        for a in range(q):
            for b in range(q):
                for c in range(q):
                    assert F.mul(a,F.mul(b,c))==F.mul(F.mul(a,b),c)
                    assert F.mul(a,F.add(b,c))==F.add(F.mul(a,b),F.mul(a,c))
        print("GF(3^%d) ok"%k)
