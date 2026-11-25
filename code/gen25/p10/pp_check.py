# G25-P10: check of Proposition PP (rank-perturbed pencil) on the exact (12,3)@6 point of G24g, plus random tests.
import sympy as sp, itertools, random
def M(entries,N=6):
    A=sp.zeros(N,N)
    for s,i,j in entries: A[i-1,j-1]=s
    return A
E=[M([(-1,2,3),(1,3,5),(1,4,1),(-1,5,2)]), M([(-1,1,5),(-1,4,4),(1,5,6),(-1,6,3)]), M([(-1,1,1),(1,2,6),(1,3,4),(1,6,2)])]
C=[M([(1,1,4),(-1,2,5),(-1,3,2),(1,5,3)]), M([(-1,3,6),(-1,4,4),(-1,5,1),(1,6,5)]), M([(-1,1,1),(1,2,6),(1,4,3),(1,6,2)])]
N=6
def swap(N):
    S=sp.zeros(N*N,N*N)
    for x in range(N):
        for y in range(N):
            S[x*N+y, y*N+x]=1   # row (x,y), col (u,v) with x=v,y=u
    return S
S=swap(N)
K=sum((sp.kronecker_product(E[l],C[l]) for l in range(3)), sp.zeros(N*N,N*N))
print("rank(Swap-K) =", (S-K).rank())
x,y,z=sp.symbols('x y z')
for name,F in (("E",E),("C",C)):
    P=x*F[0]+y*F[1]+z*F[2]
    d=sp.factor(P.det()); print(name,"det of pencil:",d)
    # min rank over the projective plane: check rank on points where minors vanish (try coordinate points and det-zero lines)
    mr=min(( (a*F[0]+b*F[1]+c*F[2]).rank() for a,b,c in itertools.product(range(-2,3),repeat=3) if (a,b,c)!=(0,0,0)))
    print(name,"min rank over small integer combinations:",mr)
