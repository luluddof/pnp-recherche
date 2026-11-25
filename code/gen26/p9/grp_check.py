# Theorem G independent tests: p >= (N^2 - sum d m^2)/2 when all G_l, G'_l commute with alpha(Gamma).
import numpy as np, itertools
from common import swap
def rk(M): return np.linalg.matrix_rank(M,tol=1e-8)
def minp(N,basis,coefs,maxq):
    """enumerate C = sum_l G_l (x) G'_l, G_l, G'_l in {sum c_i basis_i : c in coefs}, q<=maxq; return min p"""
    S=swap(N); mats=[]
    for c in itertools.product(coefs,repeat=len(basis)):
        if any(c): mats.append(sum(ci*b for ci,b in zip(c,basis)))
    prods=[np.kron(a,b) for a in mats for b in mats]
    best=N*N
    for q in range(1,maxq+1):
        for comb in itertools.combinations(range(len(prods)),q):
            K=sum(prods[i] for i in comb); best=min(best,rk(S-K))
    return best,len(prods)
# 1) Z_N regular, circulant basis {P^k}, coefficients {-1,0,1}: alpha=regular, sum d m^2 = N
for N,maxq in [(3,2),(4,1)]:
    P=np.roll(np.eye(N),1,axis=1); basis=[np.linalg.matrix_power(P,k) for k in range(N)]
    b,np_=minp(N,basis,[-1,0,1],maxq); print('Z%d circulant coeffs{-1,0,1} q<=%d: min p=%d, bound=%.1f (#prods %d)'%(N,maxq,b,(N*N-N)/2,np_))
# 2) S3 acting on R^3 by permutations = triv + std (d=1,2; m=1,1): commutant span{I, J}; bound (9-3)/2=3
J=np.ones((3,3)); basis=[np.eye(3),J]
S3_=swap(3); bb=[np.kron(a,b) for a in basis for b in basis]; vals=[-2,-1.5,-1,-2/3,-0.5,-1/3,0,1/3,0.5,2/3,1,1.5,2]
b=min(rk(S3_-sum(c*m for c,m in zip(cc,bb))) for cc in itertools.product(vals,repeat=4))
print('S3 perm rep (triv+std), K in span{I,J}^(x2) (any q) grid 13^4: min p=%d, bound=3'%b)
# 3) Z2 x Z2 ... alpha=diag(1,1,-1,-1) on R^4: commutant M2+M2 (d=1,m=2,2) bound (16-8)/2=4; random search + structured optimum
S=swap(4); rng=np.random.default_rng(0); best=16
def blk(a,b): M=np.zeros((4,4)); M[:2,:2]=a; M[2:,2:]=b; return M
for t in range(20000):
    q=rng.integers(1,4); K=0
    for l in range(q):
        A=blk(rng.integers(-1,2,(2,2)),rng.integers(-1,2,(2,2))); B=blk(rng.integers(-1,2,(2,2)),rng.integers(-1,2,(2,2)))
        K=K+np.kron(A,B)
    best=min(best,rk(S-K))
print('Z2 on R^4 (m=2,2): random q<=3 int entries: min p=%d, bound 4'%best)
# full-span optimum: K = projection of Swap onto commutant(x)commutant = sum over blocks; compute min rank achievable = theory
# 4) Z4 / Z5 / S3-regular random signed group elements, q<=4
def regS3():
    import itertools
    G=list(itertools.permutations(range(3))); idx={g:i for i,g in enumerate(G)}
    mul=lambda a,b: tuple(a[b[i]] for i in range(3))
    L=[]
    for g in G:
        M=np.zeros((6,6))
        for h in G: M[idx[mul(g,h)],idx[h]]=1
        L.append(M)
    return L, 6+4*2   # sum d^3 = 1+1+8=10
for name,(mats,sdm) in [('Z4',([np.linalg.matrix_power(np.roll(np.eye(4),1,axis=1),k) for k in range(4)],4)),
                        ('Z5',([np.linalg.matrix_power(np.roll(np.eye(5),1,axis=1),k) for k in range(5)],5)),('S3reg',(regS3()[0],10))]:
    N=mats[0].shape[0]; S=swap(N); best=N*N
    for t in range(20000):
        q=rng.integers(1,5); K=0
        for l in range(q):
            A=sum(rng.integers(-1,2)*m for m in mats) if rng.random()<0.5 else rng.choice([-1,1])*mats[rng.integers(N)]
            B=sum(rng.integers(-1,2)*m for m in mats) if rng.random()<0.5 else rng.choice([-1,1])*mats[rng.integers(N)]
            K=K+np.kron(A,B)
        best=min(best,rk(S-K))
    print(name,'random group-algebra elems q<=4: min p=%d bound=%.1f'%(best,(N*N-sdm)/2))
