# Check: Smith(M(E)) == disjoint union over k of Smith(Per_k(E))  (P6's Smith(M) lists from m_profile_out.txt, zeta=1, t=2^i)
import sys, itertools, re, ast
sys.path.insert(0,'.')
from percomp import perm, smith_vals, MOD2
from core_explore import Ecore
prof={}
for line in open('/home/user/pnp-recherche/code/gen26/p6/m_profile_out.txt'):
    m=re.match(r"a=(\d+) i=(\d+) n=\d+ Smith\(M\)=(\[.*\])",line)
    if m: prof[(int(m.group(1)),int(m.group(2)))]=ast.literal_eval(m.group(3))
for (a,i),sM in sorted(prof.items()):
    E,_=Ecore(a,1<<i,1); n=len(E)
    U=[0]
    for k in range(1,n+1):
        subs=list(itertools.combinations(range(n),k))
        P=[[perm([[E[x][y] for y in Us] for x in S]) for Us in subs] for S in subs]
        U+=smith_vals(P,600)
    U.sort()
    print(f"a={a} i={i}: Smith(M) == union_k Smith(Per_k): {U==sM}" + ("" if U==sM else f"\n  M={sM}\n  U={U}"),flush=True)
