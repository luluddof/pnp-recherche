import numpy as np, sys
N,p,q=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]); f=sys.argv[4]
z=np.load(f); n=len(z)//2; w=z[:n]+1j*z[n:]; s=N*N
G=w[:q*s].reshape(q,N,N); H=w[q*s:].reshape(q,N,N)
K=np.einsum('kxu,kyv->xyuv',G,H).reshape(s,s)
S=np.zeros((N,)*4)
for x in range(N):
  for y in range(N): S[x,y,y,x]=1
S=S.reshape(s,s)
sv=np.linalg.svd(S-K,compute_uv=False); print('sv S-K:',np.round(sv[:p+3],8))
print('ranks G:',[np.linalg.matrix_rank(g,1e-7) for g in G],'H:',[np.linalg.matrix_rank(h,1e-7) for h in H])
# pencil invariants: span{G} -- dims of common kernels/images
GG=np.hstack(G); print('rank [G1 G2 G3] (common image dim):',np.linalg.matrix_rank(GG,1e-7),' rank [G;..] (N - common kernel):',np.linalg.matrix_rank(np.vstack(G),1e-7))
HH=np.hstack(H); print('rank [H..]:',np.linalg.matrix_rank(HH,1e-7),' rank [H;..]:',np.linalg.matrix_rank(np.vstack(H),1e-7))
# eigenvalues of Theta = K S
Th=K@S; ev=np.linalg.eigvals(Th); print('eig Theta near 1:',np.sum(abs(ev-1)<1e-6),' near -1:',np.sum(abs(ev+1)<1e-6),' near 0:',np.sum(abs(ev)<1e-6))
print('other eigs:',np.round(sorted([e for e in ev if abs(e-1)>1e-6 and abs(e)>1e-6],key=abs),4))
