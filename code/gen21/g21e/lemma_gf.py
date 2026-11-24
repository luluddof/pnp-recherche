# key lemma in char 2 and 3 with G,G' drawn from a large extension field (approximates generic rank over the closure)
import numpy as np, galois, sys, itertools
q = int(sys.argv[1]); GF = galois.GF(q)
rng = np.random.default_rng(7)
rk = lambda M: int(np.linalg.matrix_rank(M))
def R(m,n): return GF.Random((m,n), seed=int(rng.integers(1<<30)))
def gen_rank(A,B,N,t=3):
    return max(rk(A@np.kron(R(N,N),R(N,N))@B) for _ in range(t))
p = GF.characteristic
P = galois.GF(p)
worst=9
for N in [2,3,4]:
    D=N*N; wN=9; wc=None
    cases=[]
    # coordinate subsets (all for N=2, random otherwise), entries in prime field
    E=np.eye(D,dtype=int)
    if N==2:
        subs=[s for r in range(1,D+1) for s in itertools.combinations(range(D),r)]
        for Rr in subs:
            for Cc in subs: cases.append((E[list(Rr)],E[:,list(Cc)]))
    else:
        for _ in range(300):
            Rr=rng.permutation(D)[:rng.integers(1,D+1)]; Cc=rng.permutation(D)[:rng.integers(1,D+1)]
            cases.append((E[Rr],E[:,Cc]))
    # random prime-field matrices of prescribed rank
    for _ in range(200 if N<4 else 60):
        a=rng.integers(1,D+1); b=rng.integers(1,D+1)
        A=(rng.integers(0,p,(a,a))@rng.integers(0,p,(a,D)))%p; B=(rng.integers(0,p,(D,b))@rng.integers(0,p,(b,b)))%p
        cases.append((A,B))
    # sym/alt spanning sets (char 2: alt = sym-diag-free)
    sym=[];alt=[]
    for i in range(N):
        for j in range(i,N):
            v=np.zeros(D,dtype=int); v[i*N+j]+=1; v[j*N+i]+=1; sym.append(v%p)
            if i<j:
                w=np.zeros(D,dtype=int); w[i*N+j]=1; w[j*N+i]=p-1; alt.append(w)
    for X in [np.array(sym),np.array(alt)]:
        for Y in [np.array(sym),np.array(alt)]: cases.append((X,Y.T))
    for A,B in cases:
        A=GF(np.asarray(A)%p); B=GF(np.asarray(B)%p)
        a=rk(A); b=rk(B)
        if a*b==0: continue
        s=gen_rank(A,B,N)
        r=s*D/(a*b)
        if r<wN: wN=r; wc=(a,b,s)
    print('GF(%d) N=%d cases=%d min s N^2/(ab)=%.3f at (a,b,s)=%s'%(q,N,len(cases),wN,wc),flush=True)
