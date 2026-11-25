import json,numpy as np, flint, itertools
from fractions import Fraction
def load(f,kH):
    d=json.load(open(f)); N=d['N']
    G=np.array([[[float(Fraction(str(v))) for v in r] for r in M] for M in d['G']]); H=np.array([[[float(Fraction(str(v))) for v in r] for r in M] for M in d[kH]])
    S=np.zeros((N,)*4)
    for x in range(N):
        for y in range(N): S[x,y,y,x]=1
    Sm=S.reshape(N*N,N*N); K=np.einsum('kxu,kyv->xyuv',G,H).reshape(N*N,N*N)
    return N,len(G),np.rint(K@Sm).astype(int)
D='/home/user/pnp-recherche/code/'
pts={'(4,4)@4':load(D+'gen20/g20b/exact_m_4_4_4.json','Gp'),'(5,3)@4':load(D+'gen20/g20b/exact_m_4_5_3.json','Gp')}
def rk(M,pr=1000003): return flint.nmod_mat([[int(v)%pr for v in r] for r in M],pr).rank()
for (a,(N1,q1,T1)),(b,(N2,q2,T2)) in itertools.combinations_with_replacement(pts.items(),2):
    T=np.kron(T1,T2).reshape(N1,N1,N2,N2,N1,N1,N2,N2).transpose(0,2,1,3,4,6,5,7).reshape((N1*N2)**2,(N1*N2)**2)
    print(a,'x',b,'-> N',N1*N2,'q',q1*q2,'p',rk(np.eye(T.shape[0],dtype=int)-T), 'p1 own',rk(np.eye(N1*N1,dtype=int)-T1))
