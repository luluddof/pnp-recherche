# Lemma: for all E,C in Mat_N: rank(Swap - E (x) C) >= N(N-1)/2.  Exact integer tests (random, low rank, nilpotent, +-I, perms).
import flint, random, itertools
def kron(A,B,N): return [[A[i][k]*B[j][l] for k in range(N) for l in range(N)] for i in range(N) for j in range(N)]
def swap(N): return [[1 if (k==j and l==i) else 0 for k in range(N) for l in range(N)] for i in range(N) for j in range(N)]
def rk(M): return flint.fmpz_mat(M).rank()
random.seed(1); worst={}
for N in range(2,7):
    S=swap(N); w=10**9
    gens=[]
    for t in range(300):
        kind=t%6
        if kind==0: E=[[random.randint(-2,2) for _ in range(N)] for _ in range(N)]; C=[[random.randint(-2,2) for _ in range(N)] for _ in range(N)]
        elif kind==1: s=random.choice([1,-1]); E=[[s*(i==j) for j in range(N)] for i in range(N)]; C=[[(i==j) for j in range(N)] for i in range(N)]
        elif kind==2:
            p=list(range(N)); random.shuffle(p); E=[[1 if p[i]==j else 0 for j in range(N)] for i in range(N)]
            C=[[1 if p[j]==i else 0 for j in range(N)] for i in range(N)]
        elif kind==3: d=[random.choice([1,-1,0,2]) for _ in range(N)]; E=[[d[i]*(i==j) for j in range(N)] for i in range(N)]; C=[[(1 if d[i]==0 else 1)*(i==j)*random.choice([1,-1]) for j in range(N)] for i in range(N)]
        elif kind==4: E=[[1 if j==i+1 else 0 for j in range(N)] for i in range(N)]; C=[[random.choice([0,1,-1]) for j in range(N)] for i in range(N)]
        else:
            u=[random.randint(-1,1) for _ in range(N)]; v=[random.randint(-1,1) for _ in range(N)]
            E=[[ (i==j)+u[i]*v[j] for j in range(N)] for i in range(N)]; C=[[ (i==j)*(-1 if i%2 else 1) for j in range(N)] for i in range(N)]
        K=kron(E,C,N); r=rk([[S[a][b]-K[a][b] for b in range(N*N)] for a in range(N*N)])
        w=min(w,r)
    print("N=%d  min rank(Swap - E(x)C) over tests = %d   N(N-1)/2 = %d  %s"%(N,w,N*(N-1)//2,"OK" if w>=N*(N-1)//2 else "VIOLATION"))
