# G27-T2 library: load exact points, build K, Theta, exact ranks.  Convention K[(x,y),(u,v)] = sum_l E_l[x,u] C_l[y,v],
# Swap[(x,y),(u,v)] = [x=v][y=u]; p = rank(Swap - K) = rank(1 - Theta), Theta = K Swap.
import json, ast, glob, numpy as np, flint
B='/home/user/pnp-recherche/code/'
SRC=[B+'gen19/g19b/exact_4_5_3.json',B+'gen19/g19b/exact_5_7_4.json',B+'gen19/g19b/exact_5_6_6.json',B+'gen19/g19b/exact_6_7_7.json',
     B+'gen24/g24g/exact_6_12_3.json']+sorted(glob.glob(B+'gen26/p11/exact_*.json'))
def load(f):
    d=json.load(open(f))
    if 'G' in d: A=[np.array([[int(x) for x in r] for r in m]) for m in d['G']]; Bm=[np.array([[int(x) for x in r] for r in m]) for m in d['Gp']]
    else: A=[np.array(m,dtype=int) for m in d['E']]; Bm=[np.array(m,dtype=int) for m in d['C']]
    return d['N'],A,Bm
def q8():
    d=json.load(open(B+'gen23/g23a/sol_Q8_8_8_8.json')); cells=ast.literal_eval(d['cells']) if isinstance(d['cells'],str) else d['cells']
    cls={}
    for a in range(8):
        for b in range(8):
            P,k,s=cells[a][b]; cls.setdefault((P,k),[]).append((a,b,s))
    G=[];H=[]
    for k,v in cls.items():
        if k[0]!='S': continue
        A=np.zeros((8,8),dtype=int)
        for a,b,s in v: A[a,b]=s
        G.append(A); H.append(A.T.copy())
    return 8,G,H
def swapm(N):
    S=np.zeros((N*N,N*N),dtype=int)
    for x in range(N):
        for y in range(N): S[x*N+y,y*N+x]=1
    return S
def Kmat(E,C): return sum(np.kron(a,b) for a,b in zip(E,C))
def frank(M): return int(np.linalg.matrix_rank(np.asarray(M,dtype=float)))
def exrank(M): return flint.fmpz_mat([[int(v) for v in r] for r in M]).rank()
def modrank(M,pr): return flint.nmod_mat([[int(v)%pr for v in r] for r in M],pr).rank()
def osr(E,C):
    N=E[0].shape[0]; R=sum(np.outer(a.ravel(),b.ravel()) for a,b in zip(E,C)); return R
def pq(N,E,C,exact=False):
    M=swapm(N)-Kmat(E,C); R=osr(E,C)
    return (exrank(M),exrank(R)) if exact else (frank(M),frank(R))
def trivial(N):
    I=np.eye(N,dtype=int); out=[]
    out.append(('I(x)I',[I],[I]))
    U=[];V=[]
    for a in range(N):
        for b in range(N):
            e=np.zeros((N,N),dtype=int); e[a,b]=1; f=np.zeros((N,N),dtype=int); f[b,a]=1; U.append(e); V.append(f)
    out.append(('Swap',U,V))
    # A_N = Swap - Phi Phi^T: osr N(N-1)/2, rank(Swap - A) = 1.  Terms: E_ab (x) E_ba - E_aa (x) E_bb  ->  a != b terms
    # realign(A) = Swap_r - I: written as sum_{a<b} (E_ab - E_ba) (x) (E_ba - E_ab) ... use exact pairs:
    if 2<=N<=4:
        U=[];V=[]
        for a in range(N):
            for b in range(a+1,N):
                e=np.zeros((N,N),dtype=int); e[a,b]=1; e[b,a]=-1; U.append(e); V.append(-e.copy())
        out.append(('A_%d'%N,U,V))
    return out
def library():
    L=[]
    for f in SRC:
        N,E,C=load(f); L.append((f.split('/')[-1],N,E,C))
    N,E,C=q8(); L.append(('Q8',N,E,C))
    return L
def save(fn,N,E,C,extra=None):
    d={"N":N,"E":[np.asarray(e).astype(int).tolist() for e in E],"C":[np.asarray(c).astype(int).tolist() for c in C],
       "convention":"K[x,y,u,v]=sum_l E_l[x,u] C_l[y,v]; Swap[x,y,u,v]=[x=v][y=u]"}
    if extra: d.update(extra)
    json.dump(d,open(fn,'w'))
if __name__=='__main__':
    for name,N,E,C in library():
        print(name,N,len(E),pq(N,E,C))
