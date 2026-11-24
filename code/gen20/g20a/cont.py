# continuation from the exact product point (t^2, N^2, N^2): warm start, then ask for smaller cut-0 rank
import numpy as np, sys, warnings, time
warnings.filterwarnings('ignore')
from trilib import *
N=int(sys.argv[1]); ps=[int(t) for t in sys.argv[2].split(',')]; tl=float(sys.argv[3]); noise=float(sys.argv[4])
b=box3([N]*4); I=np.eye(N); F=np.einsum('ad,bc->abcd',I,I); Id=np.einsum('ac,bd->abcd',I,I); Y=Id-F
L='abcdefgh'; s1=L[0]+L[3]+L[4]+L[7]; s2=L[1]+L[2]+L[5]+L[6]
mk=lambda O1,O2: np.einsum(f'{s1},{s2}->{L}',O1,O2)
A1=mk(Id,F); A2=mk(F,Y)
pr=Prob(b,ps); rng=np.random.default_rng(0); T=[]
for j,A in [(1,A1),(2,A2)]:
    U,s,Vt=np.linalg.svd(b.flat(A,j)); r=ps[j]
    k=min(r,int((s>1e-9).sum()))
    Uj=np.zeros((U.shape[0],r)); Vj=np.zeros((Vt.shape[1],r))
    Uj[:,:k]=U[:,:k]*np.sqrt(s[:k]); Vj[:,:k]=Vt[:k].T*np.sqrt(s[:k])
    T+=[Uj,Vj]
x0=pr.pack(T)
T_,R,U,s,V,E,f0=tail_state(pr,x0); print('start tail',f0,flush=True)
for i in range(int(sys.argv[5]) if len(sys.argv)>5 else 3):
    x,f,h=lm_solve(pr,x0+noise*rng.normal(size=x0.size),maxit=100000,tlim=tl)
    np.savez("cont_N%d_%s_%d.npz"%(N,"_".join(map(str,ps)),i),x=x,f=f,ps=ps)
    print(f'N={N} p={tuple(ps)} cont start {i} tail={f:.3e} its={len(h)} |x|={np.linalg.norm(x):.1f}',flush=True)
