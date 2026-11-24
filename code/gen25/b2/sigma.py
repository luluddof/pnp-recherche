# check: sigma X(x,y,u,v) -> X(v,u,y,x) fixes F, preserves Fl1-rank and Fl2-rank (osr); sigma-invariant => Pi*Fl1 symmetric
import numpy as np; rng=np.random.default_rng(1); N=4
Fl1=lambda X:X.reshape(N*N,N*N); Fl2=lambda X:X.transpose(0,2,1,3).reshape(N*N,N*N)
sig=lambda X:X.transpose(3,2,1,0)
F=np.zeros((N,)*4)
for x in range(N):
  for y in range(N): F[x,y,y,x]=1
print("sigma F == F:", np.allclose(sig(F),F))
a=rng.standard_normal((2,N*N)); b=rng.standard_normal((2,N*N))
R=(a.T@b).reshape(N,N,N,N)                       # rank 2 in Fl1
G=rng.standard_normal((3,N,N)); H=rng.standard_normal((3,N,N))
C=np.einsum('lxu,lyv->xyuv',G,H)                 # osr 3
r=np.linalg.matrix_rank
print("Fl1 rank R, sigma R:", r(Fl1(R)), r(Fl1(sig(R))), " osr C, sigma C:", r(Fl2(C)), r(Fl2(sig(C))))
Pi=Fl1(F); S=Pi@Fl1(R+sig(R)); print("Pi*Fl1(R+sigma R) symmetric:", np.allclose(S,S.T))
T=Fl2(F)@Fl2(C+sig(C)); print("Pi'*Fl2(C+sigma C) symmetric:", np.allclose(T,T.T), " Fl2(F) perm:", np.allclose(Fl2(F)@Fl2(F),np.eye(N*N)))
