# Commutative class, d=1, generic tau: which bad-edge profiles (b_1..b_l) are feasible?
# v(z) - v(z^k) = (-1)^{|z|} * 2*beta_k on good k-edges (G8d sec. 2).  beta_k = random big integers (generic).
import itertools, random, sys
from collections import Counter
l=int(sys.argv[1]); maxprof=[int(c) for c in sys.argv[2].split(',')]
n=1<<l; random.seed(1); beta=[random.randrange(10**6,10**9) for _ in range(l)]
import os
if os.environ.get("SAMETAU"): beta=[21199,12470,7310,4930][:l]   # 1/tau for tau=1,1.7,2.9,4.3 times 21199
par=[bin(z).count('1')&1 for z in range(n)]
edges=[[z for z in range(n) if not (z>>k)&1] for k in range(l)]
def consistent(good):
    # good: list of (z,k)
    adj=[[] for _ in range(n)]
    for z,k in good:
        w=z^(1<<k); g=(1-2*par[z])*2*beta[k]   # v(z)-v(w)=g
        adj[z].append((w,-g)); adj[w].append((z,g))
    v=[None]*n
    for s in range(n):
        if v[s] is not None: continue
        v[s]=0; st=[s]
        while st:
            a=st.pop()
            for b,g in adj[a]:
                if v[b] is None: v[b]=v[a]+g; st.append(b)
                elif v[b]!=v[a]+g: return False
    return True
res=Counter()
def rec(k,good,prof):
    if k==l:
        print("FOUND",prof,[e for e in good]); sys.stdout.flush(); raise SystemExit
    for nb in [maxprof[k]]:
        for B in itertools.combinations(edges[k],nb):
            Bs=set(B); g2=good+[(z,k) for z in edges[k] if z not in Bs]
            if consistent(g2): rec(k+1,g2,prof+[nb])
rec(0,[],[])
for p in sorted(res,key=lambda p:(sum(p),p)): print(p,sum(p),res[p])
