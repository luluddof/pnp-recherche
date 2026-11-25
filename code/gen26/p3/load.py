import json, numpy as np, sympy
B='/home/user/pnp-recherche/code/'
FILES=[B+'gen19/g19b/exact_4_5_3.json',B+'gen19/g19b/exact_5_7_4.json',B+'gen19/g19b/exact_5_6_6.json',
       B+'gen19/g19b/exact_6_7_7.json',B+'gen24/g24g/exact_6_12_3.json']
def load(f):
    d=json.load(open(f))
    if 'G' in d: A=[np.array(g,dtype=int) for g in [[[int(x) for x in r] for r in m] for m in d['G']]]; Bm=[np.array([[int(x) for x in r] for r in m]) for m in d['Gp']]
    else: A=[np.array(m,dtype=int) for m in d['E']]; Bm=[np.array(m,dtype=int) for m in d['C']]
    return d['N'],A,Bm
def mats(N,A,Bm):
    # rows (x,y), cols (u,v): K = sum A[x,u] B[y,v] = kron(A,B); Swap: [x=v][y=u]
    K=sum(np.kron(a,b) for a,b in zip(A,Bm))
    S=np.zeros((N*N,N*N),dtype=int)
    for x in range(N):
        for y in range(N): S[x*N+y,y*N+x]=1
    return S,K
def ismono(m): return all(np.sum(m[i]!=0)<=1 for i in range(len(m))) and all(np.sum(m[:,j]!=0)<=1 for j in range(len(m))) and set(np.unique(m))<= {-1,0,1}
def rank(M): return sympy.Matrix(M.tolist()).rank()
if __name__=='__main__':
    for f in FILES:
        N,A,Bm=load(f); S,K=mats(N,A,Bm)
        print(f.split('/')[-1],'N',N,'q',len(A),'rank(Swap-K)=',rank(S-K),'monoA',all(map(ismono,A)),'monoB',all(map(ismono,Bm)))
