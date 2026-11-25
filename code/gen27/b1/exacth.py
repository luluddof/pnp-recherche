# exact check of a rounded hada solution: entries -> rationals (limit_denominator), rebuild K in original basis, flint ranks
import numpy as np, sys, flint
from fractions import Fraction
from scipy.linalg import hadamard
fn=sys.argv[1]; N=int(sys.argv[2]); md=int(sys.argv[3]) if len(sys.argv)>3 else 1000
z=np.load(fn); fr=[Fraction(float(v)).limit_denominator(md) for v in z]
print('max rounding err',max(abs(float(f)-v) for f,v in zip(fr,z)))
X=np.array(fr[:N*N],dtype=object).reshape(N,N); Y=np.array(fr[N*N:],dtype=object).reshape(N,N)
Hm=hadamard(N).astype(object)
A=Hm.dot(X).dot(Hm)/(N*N); B=Hm.dot(Y).dot(Hm)/N
s=N*N; K=[[Fraction(0)]*s for _ in range(s)]; KR=[[Fraction(0)]*s for _ in range(s)]
for g in range(N):
    for x in range(N):
        for u in range(N):
            a=A[x^g,u^g]
            if a==0: continue
            for y in range(N):
                for v in range(N):
                    b=B[y^g,v^g]
                    if b==0: continue
                    K[x*N+y][u*N+v]+=a*b; KR[x*N+u][y*N+v]+=a*b
for x in range(N):
    for y in range(N): K[x*N+y][y*N+x]-=1   # K - S
den=1
fq=lambda M: flint.fmpq_mat([[flint.fmpq(e.numerator,e.denominator) for e in r] for r in M])
print('exact rank_Q(S-K) =',fq(K).rank(),' exact osr_Q(K) =',fq(KR).rank(),flush=True)
for pr in [1000003, 998244353]:
    def m(e): return (e.numerator%pr)*pow(e.denominator%pr,-1,pr)%pr
    print(pr, flint.nmod_mat([[m(e) for e in r] for r in K],pr).rank(), flint.nmod_mat([[m(e) for e in r] for r in KR],pr).rank())
