# Leaf S2.1.1b: monomial-C model. K = sum_{l<q} E_l (x) C_l, E_l, C_l signed partial permutation matrices,
# R = Swap - K arbitrary; p = rank(Swap - K). Simulated annealing on F = sum of singular values beyond index ptarget.
# convention K[(x,y),(u,v)] = E[x,u] C[y,v], Swap[(x,y),(u,v)] = [x=v][y=u]
import numpy as np, sys, time
N=int(sys.argv[1]); q=int(sys.argv[2]); pt=int(sys.argv[3]); secs=float(sys.argv[4]); seed=int(sys.argv[5])
allowzero=1 if len(sys.argv)<7 else int(sys.argv[6])
rng=np.random.default_rng(seed)
n=N*N
S=np.zeros((n,n))
for x in range(N):
  for y in range(N): S[x*N+y,y*N+x]=1
def mat(pi,s):
  M=np.zeros((N,N)); M[np.arange(N),pi]=s; return M
def K(st):
  return sum(np.kron(mat(*st[2*l]),mat(*st[2*l+1])) for l in range(q))
def F(st):
  sv=np.linalg.svd(S-K(st),compute_uv=False)
  return sv[pt:].sum(), int((sv>1e-8).sum())
def rand():
  return [(rng.permutation(N), rng.choice([-1,1],N).astype(float)) for _ in range(2*q)]
def move(st):
  st=[(a.copy(),b.copy()) for a,b in st]
  i=rng.integers(2*q); pi,s=st[i]
  if rng.random()<0.5:
    a,b=rng.choice(N,2,replace=False); pi[a],pi[b]=pi[b],pi[a]
  else:
    a=rng.integers(N); choices=[-1.,1.]+([0.] if allowzero else []); s[a]=rng.choice([c for c in choices if c!=s[a]])
  return st
t0=time.time(); best=(1e9,None,None)
while time.time()-t0<secs:
  st=rand(); f,rk=F(st); T=1.0; it=0
  while time.time()-t0<secs and it<20000:
    it+=1; st2=move(st); f2,rk2=F(st2)
    if f2<=f or rng.random()<np.exp((f-f2)/T): st,f,rk=st2,f2,rk2
    T*=0.9995
    if f<best[0]: best=(f,rk,st)
    if f<1e-9: break
  if best[0]<1e-9: break
f,rk,st=best
print("N",N,"q",q,"ptarget",pt,"bestF",f,"rank",rk,"time",round(time.time()-t0))
if f<1e-9:
  # exact rank check over Q via sympy-free integer elimination mod two primes
  M=(S-K(st)).astype(np.int64)
  def rk_mod(M,P):
    A=M%P; r=0; A=A.copy(); rows,cols=A.shape
    for c in range(cols):
      piv=None
      for i in range(r,rows):
        if A[i,c]: piv=i;break
      if piv is None: continue
      A[[r,piv]]=A[[piv,r]]; inv=pow(int(A[r,c]),P-2,P); A[r]=(A[r]*inv)%P
      for i in range(rows):
        if i!=r and A[i,c]: A[i]=(A[i]-A[i,c]*A[r])%P
      r+=1
    return r
  print("rank mod 10007:",rk_mod(M,10007),"mod 3:",rk_mod(M,3), "mod 2:", rk_mod(M,2))
  print([ (list(a),list(b)) for a,b in st])
