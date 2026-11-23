# Check: OpenAI local deformation J(t)=p(sum (1-t/j)X_j) + E, E = X^{(x)l} - p(sum X_j).
# Lemma (G4c): any family with deg_t J <= d, J(0)=X^{(x)l}, J(i)=K_i+E (K_i preserves bit i, E common)
# forces rank(E) >= ball(l, floor((l-d-1)/2)). Compare with rank of OpenAI's E.
import numpy as np, itertools
from fractions import Fraction
from math import comb
def ball(l,r): return sum(comb(l,j) for j in range(0,r+1)) if r>=0 else 0
def openai_E_rank(l,d):
    # window of d+1 consecutive k centred at l/2; p(l-2k)=(-1)^k there
    c=l/2; ks=sorted(range(l+1),key=lambda k:(abs(k-c),k))[:d+1]
    pts=[(l-2*k,(-1)**k) for k in ks]
    def p(z):
        s=Fraction(0)
        for i,(zi,yi) in enumerate(pts):
            t=Fraction(yi)
            for j,(zj,_) in enumerate(pts):
                if j!=i: t*=Fraction(z-zj,zi-zj)
            s+=t
        return s
    # rank(E) = #{x in {-1,1}^l : p(sum x) != prod x}  (Hadamard diagonalisation)
    return sum(comb(l,k) for k in range(l+1) if p(l-2*k)!=(-1)**k)
def brute(l,d):
    # build actual matrices and check rank numerically + niceness at t=i
    n=2**l; X=np.array([[0,1],[1,0]]); I=np.eye(2)
    Xs=[]
    for j in range(l):
        M=np.array([[1.]])
        for k in range(l): M=np.kron(M, X if k==j else I)
        Xs.append(M)
    Xl=np.array([[1.]])
    for k in range(l): Xl=np.kron(Xl,X)
    c=l/2; ks=sorted(range(l+1),key=lambda k:(abs(k-c),k))[:d+1]
    zs=[l-2*k for k in ks]; ys=[(-1)**k for k in ks]
    coef=np.polyfit(zs,ys,d)
    def pm(A):
        R=np.zeros_like(A); P=np.eye(n)
        for a in coef[::-1]: R+=a*P; P=P@A
        return R
    E=Xl-pm(sum(Xs))
    return np.linalg.matrix_rank(E,tol=1e-6)
print(" l  d  OpenAI_rank(E)  lower_bound  ratio  [brute rank for l<=8]")
for l in [6,8,10,16,24,40,64,100]:
    for d in range(0,l,max(1,l//6)):
        ub=openai_E_rank(l,d); lb=ball(l,(l-d-1)//2)
        br=brute(l,d) if l<=8 else ''
        print(f"{l:3d} {d:3d} {ub:14d} {lb:12d} {ub/max(lb,1):6.2f}  {br}")
