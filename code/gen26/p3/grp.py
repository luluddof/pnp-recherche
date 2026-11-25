# Theorem G check: A_l=+-L(g_l), B_l=+-L(h_l) left-regular in a group; claim p >= (N^2 - sum_rho d_rho^3)/2
import numpy as np, itertools, random, sys
from sympy.combinatorics import Permutation, PermutationGroup
from sympy.combinatorics.named_groups import SymmetricGroup, DihedralGroup, CyclicGroup, AlternatingGroup
def quat():
    # Q8 as permutations of 8 elements via left mult; build via sympy from generators i,j on regular rep
    # elements: +-1,+-i,+-j,+-k encoded 0..7 ; multiplication table
    names=['1','i','j','k']; tab={('1',x):(1,x) for x in names}
    for x in names: tab[(x,'1')]=(1,x)
    for x in 'ijk': tab[(x,x)]=(-1,'1')
    tab[('i','j')]=(1,'k');tab[('j','i')]=(-1,'k');tab[('j','k')]=(1,'i');tab[('k','j')]=(-1,'i');tab[('k','i')]=(1,'j');tab[('i','k')]=(-1,'j')
    els=[(s,x) for s in (1,-1) for x in names]; idx={e:t for t,e in enumerate(els)}
    def mul(a,b): s,x=tab[(a[1],b[1])]; return (a[0]*b[0]*s,x)
    return [[idx[mul(a,b)] for b in els] for a in els]   # table[a][b]=a*b
def table_from_group(G):
    els=list(G.elements); idx={e:t for t,e in enumerate(els)}
    return [[idx[a*b] for b in els] for a in els]
def irrep_dims(tab):
    # class sizes -> use character-free method: dims from decomposition of group algebra center? use numeric: sum d^3 via
    # dims of irreps computed from eigen decomposition of regular rep commutant is heavy; instead brute: dimension of
    # algebra of class sums -> number of irreps k; we compute d's by numerically block-diagonalising a random element of C[G] x C[G]^op
    n=len(tab); L=[np.zeros((n,n)) for _ in range(n)]
    for a in range(n):
        for b in range(n): L[a][tab[a][b],b]=1
    R=[np.zeros((n,n)) for _ in range(n)]
    for a in range(n):
        for b in range(n): R[a][tab[b][a],b]=1
    # random element of center? use: the commutant of L is R-algebra; isotypic decomposition via random symmetric element of
    # span{L_g + L_g^T} combined with R : eigenvalues of X = sum c_g L_g (sym) restricted... simpler: d_rho^2 multiplicities of
    # eigenvalues of random central element Z = sum_C c_C sum_{g in C} L_g  (class function, real-symmetrised)
    inv=[ [b for b in range(n) if tab[a][b]==0][0] for a in range(n)]
    classes=[];seen=set()
    for a in range(n):
        if a in seen: continue
        cl=set(tab[tab[g][a]][inv[g]] for g in range(n)); seen|=cl; classes.append(sorted(cl))
    c=np.random.randn(len(classes)); Z=sum(c[i]*sum(L[g] for g in cl) for i,cl in enumerate(classes))
    Z=Z+np.array([[0]]) ; ev=np.linalg.eigvals(Z)
    # group eigenvalues (central element acts as scalar on each isotypic block of size d^2)
    ev=sorted(ev,key=lambda z:(round(z.real,6),round(z.imag,6))); groups=[]
    for z in ev:
        for gr in groups:
            if abs(gr[0]-z)<1e-6: gr.append(z); break
        else: groups.append([z])
    ds=[int(round(np.sqrt(len(g)))) for g in groups]
    return ds,L
def test(tab,q,trials,name):
    ds,L=irrep_dims(tab); n=len(tab); N=n
    bound=(N*N-sum(d**3 for d in ds))/2
    S=np.zeros((N*N,N*N))
    for x in range(N):
        for y in range(N): S[x*N+y,y*N+x]=1
    best=10**9
    for t in range(trials):
        K=sum(random.choice([1,-1])*np.kron(L[random.randrange(n)],L[random.randrange(n)]) for _ in range(q))
        p=np.linalg.matrix_rank(S-K); best=min(best,p)
        assert p>=bound-1e-9, (name,p,bound)
    print(name,'N',N,'irrep dims',sorted(ds),'q',q,'bound (N^2-sum d^3)/2 =',bound,'min p over',trials,'random =',best)
random.seed(1)
G={'Z6':table_from_group(CyclicGroup(6)),'S3':table_from_group(SymmetricGroup(3)),'Q8':quat(),'D4':table_from_group(DihedralGroup(4)),
   'Z2^3':table_from_group(PermutationGroup([Permutation([1,0,2,3,4,5]),Permutation([0,1,3,2,4,5]),Permutation([0,1,2,3,5,4])])),
   'A4':table_from_group(AlternatingGroup(4)),'D6':table_from_group(DihedralGroup(6))}
if __name__=="__main__":
  for name,tab in G.items():
    for q in (1,2,3,6):
        test(tab,q,60,name)
