# frontier table N = 2..12 from all exact points (repo + cons/ + exact/), mirror symmetry (p,q) <-> (q,p), trend fits.
import glob, os, json, math, numpy as np
from lib import SRC
pts={}
def addp(N,p,q,src):
    pts.setdefault(N,[]).append((p,q,src))
for f in SRC+glob.glob('cons/c_*.json')+glob.glob('exact/exact_*.json'):
    b=os.path.basename(f).replace('.json','').split('_'); N,p,q=int(b[-3]),int(b[-2]),int(b[-1]); addp(N,p,q,f)
# analytic points: f(N,1) = N(N-1)/2, (1, N(N-1)/2) mirror; Q8 (8,8)@8, (8,8)@7 via deletion (in cons)
for N in range(2,13): addp(N,N*(N-1)//2,1,'I(x)I')
def frontier(N):
    P=set((p,q) for p,q,_ in pts.get(N,[]) if p>0 and q>0)
    P|=set((q,p) for p,q in P)           # mirror (realignment) symmetry
    # monotonicity in N: restriction f(N,q) <= f(N+1,q)
    for M in range(N+1,13):
        P|=set((p,q) for p,q,_ in pts.get(M,[]) if p>0 and q>0)|set((q,p) for p,q,_ in pts.get(M,[]) if p>0 and q>0)
    Fr=[]
    for p,q in sorted(P,key=lambda t:(t[1],t[0])):
        if not any((p2<=p and q2<=q and (p2,q2)!=(p,q)) for p2,q2 in P): Fr.append((p,q))
    return sorted(set(Fr),key=lambda t:t[1])
rows=[]
print("N | frontier (p,q), q<=p half | b(N) | min(th+T) all | min(th+T) p,q>=3 | min pq/N^2 (p,q>=3) | min th^2+T^2 (p,q>=3)")
for N in range(2,13):
    Fr=frontier(N); L=math.log(N)
    half=[(p,q) for p,q in Fr if q<=p]
    b=min(max(p,q) for p,q in Fr)
    s_all=min(math.log(p*q)/L for p,q in Fr)
    inner=[(p,q) for p,q in Fr if p>=3 and q>=3]
    s_in=min((math.log(p*q)/L,(p,q)) for p,q in inner) if inner else None
    c_in=min((p*q/N**2,(p,q)) for p,q in inner) if inner else None
    r2=min(((math.log(p)/L)**2+(math.log(q)/L)**2,(p,q)) for p,q in inner) if inner else None
    rows.append(dict(N=N,frontier=half,b=b,s_all=s_all,s_in=s_in,c_in=c_in,r2=r2))
    print(N,'|',' '.join(f"({p},{q})" for p,q in half),'| b =',b,f'| {s_all:.3f} |',f"{s_in[0]:.3f} {s_in[1]}" if s_in else '-', '|',f"{c_in[0]:.2f} {c_in[1]}" if c_in else '-','|',f"{r2[0]:.3f} {r2[1]}" if r2 else '-')
json.dump(rows,open('trend_table.json','w'),default=str,indent=0)
