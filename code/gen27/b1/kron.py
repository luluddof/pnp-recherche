# KT1: Kronecker products of exact optima: Theta = Theta1 (x) Theta2 (cells (x1x2, y1y2)), q = q1 q2; p = rank(1-Theta) exact mod prime.
import numpy as np, json, glob, itertools, flint
D='/home/user/pnp-recherche/code/gen26/p11/'
def load(f):
    d=json.load(open(f)); G=np.array(d['E']); H=np.array(d['C']); N=d['N']
    S=np.zeros((N,)*4,dtype=int)
    for x in range(N):
        for y in range(N): S[x,y,y,x]=1
    K=np.einsum('kxu,kyv->xyuv',G,H); Th=(K.reshape(N*N,N*N)@S.reshape(N*N,N*N))
    return N,d['rank'],d['osr'],Th,G,H
def rk(M,pr=1000003): return flint.nmod_mat([[int(v)%pr for v in r] for r in M],pr).rank()
fs=sorted(glob.glob(D+'exact_*.json'))
pts=[load(f) for f in fs if any(t in f for t in ['6_8_6','7_20_3','7_11_7','6_8_6'])]
pts=[load(D+'exact_6_8_6.json'),load(D+'exact_7_11_7.json'),load(D+'exact_7_20_3.json')]
for (N1,p1,q1,T1,G1,H1),(N2,p2,q2,T2,G2,H2) in itertools.combinations_with_replacement(pts,2):
    # Theta on cells (x1,y1,x2,y2) ordering -> reorder to ((x1,x2),(y1,y2))
    T=np.kron(T1,T2).reshape(N1,N1,N2,N2,N1,N1,N2,N2).transpose(0,2,1,3,4,6,5,7).reshape((N1*N2)**2,(N1*N2)**2)
    n=(N1*N2)**2; p=rk(np.eye(n,dtype=int)-T)
    ev1=np.round(np.linalg.eigvals(T1).real,6); ev2=np.round(np.linalg.eigvals(T2).real,6)
    N=N1*N2
    print(f"({p1},{q1})@{N1} x ({p2},{q2})@{N2}: N={N} q={q1*q2} p={p}  [p/N^2={p/n:.3f}; naive p1 n2+p2 n1-p1p2={p1*N2**2+p2*N1**2-p1*p2}; logN p={np.log(p)/np.log(N):.3f} logN q={np.log(q1*q2)/np.log(N):.3f}]",flush=True)
