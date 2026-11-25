# check Lemma MC2 (i): cells on triangles of the undirected cell graph <= 8 q^3 N; also try structured (group) perms
import numpy as np, itertools, sys
def run(N,q,sig,tau):
  nb=[set() for _ in range(N*N)]
  for l in range(q):
    for a in range(N):
      for b in range(N):
        c=a*N+b; d=sig[l][b]*N+tau[l][a]
        if d!=c: nb[c].add(d); nb[d].add(c)
  tri=set()
  for c in range(N*N):
    for d in nb[c]:
      for e in nb[c]&nb[d]: tri.update((c,d,e))
  return len(tri)
rng=np.random.default_rng(0)
for N,q in [(20,2),(20,3),(30,3),(30,5),(40,4)]:
  sig=[rng.permutation(N) for _ in range(q)]; tau=[rng.permutation(N) for _ in range(q)]
  t=run(N,q,sig,tau); print("random",N,q,"tri cells",t,"of",N*N,"bound",8*q**3*N)
# structured: cyclic shifts sigma_l(b)=b+s_l, tau_l(a)=a+t_l mod N
for N,q in [(20,3),(30,4)]:
  sig=[(np.arange(N)+rng.integers(N))%N for _ in range(q)]; tau=[(np.arange(N)+rng.integers(N))%N for _ in range(q)]
  t=run(N,q,sig,tau); print("cyclic",N,q,"tri cells",t,"of",N*N,"bound",8*q**3*N)
# adversarial: sigma_l = tau_l^{-1} (P_l involutions) with commuting shifts
N,q=30,4
sh=[0,1,2,3]; sig=[(np.arange(N)+s)%N for s in sh]; tau=[(np.arange(N)-s)%N for s in sh]
print("involutive shifts",N,q,"tri cells",run(N,q,sig,tau),"of",N*N)
