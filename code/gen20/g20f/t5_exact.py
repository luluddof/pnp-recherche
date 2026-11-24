import numpy as np, sympy as sp, random
def Smat(k,m):
    S=np.zeros((m*k,k*m),dtype=np.int64)
    for i in range(k):
        for j in range(m): S[j*k+i, i*m+j]=1
    return S
# exact power-sum check (determines charpoly in char 0)
ok=True
for (k,m) in [(1,3),(2,3),(3,2),(2,4),(3,5),(4,2),(3,3),(4,4)]:
    for trial in range(3):
        gt=sp.Matrix(m,k,lambda a,b:random.randint(-3,3)); ga=sp.Matrix(k,m,lambda a,b:random.randint(-3,3))
        Th=sp.Matrix(Smat(k,m).tolist()).inv()*sp.kronecker_product(gt,ga); P=ga*gt
        T=sp.eye(k*m)
        for n in range(1,k*m+1):
            T=T*Th
            pred=(P**((n)//2)).trace()**2 if n%2==0 else (P**n).trace()
            if T.trace()!=pred: ok=False; print('FAIL',k,m,n)
print('exact power sums tr Theta^n == formula for n<=km, all cases:',ok)
# Theorem 5 as literally stated, q=1, N=2, k=3 (k > N): bound vs feasible (1,1)
N=2; F=sp.zeros(4,4)
for x in range(2):
    for y in range(2): F[x*2+y,y*2+x]=1
print('N=2: rank(F - I(x)I) =',(F-sp.eye(4)).rank(),'(so (p,q)=(1,1) feasible)')
for k in range(1,5):
    m=N-k*0; kp=min(k,m); print(f'  k={k}: Thm5 literal bound p >= {k*m-kp*(kp+1)//2}')
