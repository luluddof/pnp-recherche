# [E] exact verification: round a sparsified solution, check rank_Q(Swap - K) and osr(K) with flint; print as matrix units.
import numpy as np, flint, sys, json
N,p,q=map(int,sys.argv[1].split(',')); fn=sys.argv[2]
z=np.load(fn); s=N*N
G=np.rint(z[:q*s].reshape(q,N,N)).astype(int); H=np.rint(z[q*s:].reshape(q,N,N)).astype(int)
print("max rounding error:", np.abs(z-np.concatenate([G.ravel(),H.ravel()])).max())
K=np.einsum('kxu,kyv->xyuv',G,H); S=np.zeros((N,)*4,dtype=int)
for x in range(N):
    for y in range(N): S[x,y,y,x]=1
M=(S-K).reshape(s,s); R=np.einsum('kxu,kyv->xuyv',G,H).reshape(s,s)
fz=lambda A: flint.fmpz_mat([[int(v) for v in r] for r in A])
rk=fz(M).rank(); osr=fz(R).rank()
print(f"N={N}: exact rank(Swap-K) = {rk}, exact osr(K) = {osr}, N(N-1)/2 = {N*(N-1)//2}")
for pr in [2,3,5,7]:
    print(f"  rank mod {pr}: {flint.nmod_mat([[int(v)%pr for v in r] for r in M],pr).rank()}, osr mod {pr}: {flint.nmod_mat([[int(v)%pr for v in r] for r in R],pr).rank()}")
show=lambda A:" ".join(f"{'+' if A[a,b]>0 else '-'}{'' if abs(A[a,b])==1 else abs(A[a,b])}E{a+1}{b+1}" for a in range(N) for b in range(N) if A[a,b]!=0)
for l in range(q): print(f"  E_{l+1} = {show(G[l])}    C_{l+1} = {show(H[l])}")
json.dump({"N":N,"rank":rk,"osr":osr,"E":G.tolist(),"C":H.tolist(),"convention":"K[x,y,u,v]=sum_l E_l[x,u] C_l[y,v]; Swap[x,y,u,v]=[x=v][y=u]"},open(f"exact_{N}_{rk}_{osr}.json","w"))
