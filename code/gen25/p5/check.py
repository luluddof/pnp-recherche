# Exact checks of the two core identities, translated into the project's cut-rank language.
import itertools, numpy as np, sympy as sp
from fractions import Fraction
from math import comb

import flint
def rank_Q(M):
    conv=lambda x: flint.fmpq(int(Fraction(x).numerator),int(Fraction(x).denominator))
    return flint.fmpq_mat([[conv(x) for x in row] for row in M]).rank()

# ---------- (1) Fischer Lemma 3.1 swap identity on [q]^2 x [q]^2 ----------
# sum_ij tr(U_ij V_ij) = sum_ij tr(U_ij V_ji) + sum_{i<j} tr((U_ij-U_ji)(V_ij-V_ji))
# As 4-index tensors T[(i,j),(k,l)] (U index (i,j), V index (k,l)):
#   Id = Swap + sum_{i<j} a_ij (x) a_ij,  a_ij = e_ij - e_ji
def fischer(q):
    n=q*q; idx=lambda i,j:i*q+j
    Id=np.zeros((n,n),int); Sw=np.zeros((n,n),int)
    for i in range(q):
        for j in range(q):
            Id[idx(i,j),idx(i,j)]=1; Sw[idx(i,j),idx(j,i)]=1
    corr=[]
    for i in range(q):
        for j in range(i+1,q):
            a=np.zeros(n,int); a[idx(i,j)]=1; a[idx(j,i)]=-1
            corr.append(np.outer(a,a))
    assert (Id==Sw+sum(corr)).all()
    # cut A = (U-index)|(V-index) = (i,j)|(k,l);  cut B = (i,k)|(j,l)
    def cutB(T):
        T4=T.reshape(q,q,q,q)            # i j k l
        return T4.transpose(0,2,1,3).reshape(n,n)
    rk=lambda T:(np.linalg.matrix_rank(T),np.linalg.matrix_rank(cutB(T)))
    return rk(Id),rk(Sw),rk(corr[0]),rk(sum(corr))
for q in (2,3,4,5):
    print("Fischer q=%d  ranks (cutA,cutB): Id %s Swap %s one-correction %s all-corrections %s  #corr=%d"%((q,)+fischer(q)+(comb(q,2),)))

# ---------- (2) CIMY Lemma 2.4 (Hurwitz/Hrubes identity) ----------
P={'I':np.eye(2,dtype=int),'X':np.array([[0,1],[1,0]]),'Z':np.array([[1,0],[0,-1]]),'J':np.array([[0,1],[-1,0]])}
def anti(s,t):  # do pauli strings anticommute
    c=sum(1 for x,y in zip(s,t) if x!='I' and y!='I' and x!=y); return c%2==1
def mat(s):
    M=np.array([[1]])
    for ch in s: M=np.kron(M,P[ch])
    return M
def family(r,m):
    cands=[''.join(s) for s in itertools.product('IXZJ',repeat=r) if s.count('J')%2==1]
    best=[]
    def rec(cur,start):
        nonlocal best
        if len(cur)>len(best): best=cur[:]
        if len(best)>=m: return True
        for i in range(start,len(cands)):
            if all(anti(cands[i],c) for c in cur):
                cur.append(cands[i])
                if rec(cur,i+1): return True
                cur.pop()
        return False
    rec([],0); return [mat(s) for s in best[:m]]
def cimy(m,k):
    r=m//2; nA=2**r; E=family(r,m); assert len(E)==m
    for a in range(m):
        assert (E[a].T==-E[a]).all() and (E[a]@E[a].T==np.eye(nA)).all()
        for b in range(a+1,m): assert (E[a]@E[b]==-E[b]@E[a]).all()
    B=list(itertools.combinations(range(m),k)); nB=len(B)
    EB=[]
    for beta in B:
        M=np.eye(nA,dtype=int)
        for i in beta: M=M@E[i]
        EB.append(M)
    # p(z) = prod_{even i<k}(z-i)/K ; c_i via inversion c_i = sum_j (-1)^{i-j} C(i,j) p(j)
    K=1
    for i in range(2,k+1,2): K*=i
    pz=lambda z: Fraction(np.prod([z-i for i in range(0,k,2)]),K)
    c=[sum((-1)**(i-j)*comb(i,j)*pz(j) for j in range(i+1)) for i in range(k//2+1)]
    Lam=[l for s in range(k//2+1) for l in itertools.combinations(range(m),s)]
    L=nA*len(Lam)
    # h^{(j,lam)}_{a,beta} = E_beta[a,j] [lam subset beta], w = c_|lam|
    D=nA*nB
    Hm=np.full((D,L),0,dtype=object); w=[]
    col=0
    for j in range(nA):
        for lam in Lam:
            w.append(Fraction(c[len(lam)]))
            for bi,beta in enumerate(B):
                if set(lam)<=set(beta):
                    for a in range(nA):
                        Hm[a*nB+bi,col]=int(EB[bi][a,j])
            col+=1
    from math import lcm
    den=1
    for x in w: den=lcm(den,Fraction(x).denominator)
    wi=np.array([int(Fraction(x)*den) for x in w],dtype=np.int64)
    Hm=Hm.astype(np.int64)
    H=(Hm*wi)@Hm.T                     # = den * H_true
    H4=H.reshape(nA,nB,nA,nB); Ht=H4.transpose(0,3,2,1).reshape(D,D)
    ok=bool(((H+Ht)==2*den*np.eye(D,dtype=np.int64)).all())
    rk=lambda A: flint.fmpz_mat(A.tolist()).rank()
    rH=rk(H); rHt=rk(Ht)
    return nA,nB,L,ok,rH,rHt
for (m,k) in [(8,2),(8,4)]:
    nA,nB,L,ok,rH,rHt=cimy(m,k)
    M=nA*nB
    print("CIMY m=%d k=%d: nA=%d nB=%d  identity(2) exact: %s  L=%d  rank_cut1(H)=%d  rank_cut2(H)=%d  full=%d  exponent log rank/log full = %.3f"%(m,k,nA,nB,ok,L,rH,rHt,M,np.log(rH)/np.log(M)))
