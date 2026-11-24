"""Exhaustive lambda=1 over GF(2): all 2^16 two-piece identities R + C = X (C = X+R).
Tensor-power pieces A_i=[T^i](R+TC)^{(x)t}; cost_t = max_i log2 rank_row(A_i)+log2 rank_col(A_i), /(4t)."""
import numpy as np, itertools, sys
from math import log2
from gf2asym import gf2rank
def pack(M):
    return [int("".join(map(str,row)),2) for row in (M%2).astype(int)]
def rk(M): return gf2rank(pack(M))
I4=np.eye(4,dtype=np.int64)
polys=[q for q in itertools.product(range(2),repeat=5) if any(q)]
def key(A):
    pw=[I4]; 
    for _ in range(4): pw.append(pw[-1]@A%2)
    return tuple(rk(sum(c*P for c,P in zip(q,pw))%2) for q in polys)
def pieces(A,t):
    Cm=(I4+A)%2
    P=[np.array([[1]],dtype=np.int64)]
    for s in range(t):
        new=[np.zeros((4**(s+1),4**(s+1)),dtype=np.int64) for _ in range(len(P)+1)]
        for i,Q in enumerate(P):
            new[i]+=np.kron(Q,A); new[i+1]+=np.kron(Q,Cm)
        P=[x%2 for x in new]
    return P
cache={}
def rv(A,t):
    k=(key(A),t)
    if k not in cache: cache[k]=[rk(Q) for Q in pieces(A,t)]
    return cache[k]
pairs={}
for bits in range(1<<16):
    R=np.array([(bits>>j)&1 for j in range(16)],dtype=np.int64).reshape(2,2,2,2)
    A=R.transpose(0,1,3,2).reshape(4,4); B=R.transpose(0,2,3,1).reshape(4,4)
    k=(key(A),key(B))
    if k not in pairs: pairs[k]=(bits,A,B)
print("distinct (class_row,class_col) pairs:",len(pairs)); sys.stdout.flush()
def cost(A,B,t):
    ra,rb=rv(A,t),rv(B,t)
    return max((log2(a)+log2(b)) if a and b else 0 for a,b in zip(ra,rb))/(4*t), ra, rb
res=[]
for k,(bits,A,B) in pairs.items():
    c,ra,rb=cost(A,B,4); res.append((c,bits,ra,rb))
res.sort()
print("t=4: best 12 of",len(res))
for c,bits,ra,rb in res[:12]: print("  rho_4=%.4f R=%04x row %s col %s"%(c,bits,ra,rb))
from collections import Counter
print("histogram of rho_4 (rounded):",sorted(Counter(round(r[0],3) for r in res).items())[:8])
print("t=5 for the best 6:")
for c,bits,ra,rb in res[:6]:
    _,A,B=[p for p in pairs.values() if p[0]==bits][0]
    c5,ra5,rb5=cost(A,B,5); print("  rho_5=%.4f R=%04x row %s col %s"%(c5,bits,ra5,rb5))
