import json, numpy as np, sys, itertools
D='/home/user/pnp-recherche/code/'
files=[D+'gen19/g19b/exact_4_5_3.json',D+'gen20/g20b/exact_m_4_5_3.json',D+'gen20/g20b/exact_m_4_4_4.json',D+'gen19/g19b/exact_5_7_4.json',D+'gen19/g19b/exact_5_6_6.json',D+'gen19/g19b/exact_6_7_7.json']
def load(f):
    d=json.load(open(f)); N=d['N']
    G=np.array([[[int(e) for e in r] for r in l] for l in d['G']]); Gp=np.array([[[int(e) for e in r] for r in l] for l in d['Gp']])
    return N,G,Gp
def Rmat(N,G,Gp):
    C=np.einsum('lxu,lyv->xyuv',G,Gp)
    X=np.zeros((N,)*4,int)
    for x in range(N):
        for y in range(N): X[x,y,y,x]=1
    return (X-C).reshape(N*N,N*N), C
if __name__=="__main__":
  for f in files:
    N,G,Gp=load(f); q=len(G)
    R,C=Rmat(N,G,Gp)
    print("=====",f.split('/')[-1],"N",N,"q",q,"rank R",np.linalg.matrix_rank(R))
    col=np.full((N,N),'.',dtype=object); colp=np.full((N,N),'.',dtype=object)
    for l in range(q):
        for x in range(N):
            for u in range(N):
                if G[l,x,u]: col[x,u]=('+' if G[l,x,u]>0 else '-')+str(l)
                if Gp[l,u,x]: colp[x,u]=('+' if Gp[l,u,x]>0 else '-')+str(l)   # G'^T
    print("G (colour/sign at [x,u]):"); 
    for r in col: print('  '+' '.join(f"{e:>3}" for e in r))
    print("G'^T:")
    for r in colp: print('  '+' '.join(f"{e:>3}" for e in r))
    same=all((G[l]==Gp[l].T).all() for l in range(q)); print("G'_l == G_l^T for all l:",same)
    # kernel of R (W with W^T = sum G W G'^T)
    u,s,vt=np.linalg.svd(R); rk=(s>1e-9).sum()
    print("singular values of R:",np.round(s[:rk],3))
    # diag coverage
    S=[(x,y) for x in range(N) for y in range(N) if C[x,y,y,x]!=0]
    print("|S| =",len(S),"C(x,y,y,x) values:",sorted(set(C[x,y,y,x] for x in range(N) for y in range(N))))
    print("nnz rows of R:",(np.abs(R).sum(1)>0).sum(),"nnz cols:",(np.abs(R).sum(0)>0).sum(), "nnz(R)",(R!=0).sum(),"nnz(C)",(C!=0).sum())
