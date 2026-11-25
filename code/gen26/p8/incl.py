# Is rowspace(C_k(E)) contained in rowspace(Per_k(E)) over O (=> C_k = A Per_k, A integral => (A))? Also column version.
import sys, itertools
sys.path.insert(0,'.')
from percomp import perm, smith_vals, MOD2
from union_check2 import Eblk
from fractions import Fraction
def det(M):
    n=len(M)
    if n==0: return 1
    tot=0
    for p in itertools.permutations(range(n)):
        s=1
        for i in range(n):
            for j in range(i+1,n):
                if p[i]>p[j]: s=-s
        pr=1
        for i in range(n): pr=pr*M[i][p[i]]%MOD2
        tot+=s*pr
    return tot%MOD2
for (a,t,z) in [(3,1,1),(3,2,1),(3,-3,Fraction(1,3)),(4,1,1),(4,2,1),(4,4,-3)]:
    E=Eblk(a,t,z,0); n=len(E)
    for k in range(2,min(n,4)+1):
        subs=list(itertools.combinations(range(n),k))
        P=[[perm([[E[x][y] for y in U] for x in S]) for U in subs] for S in subs]
        C=[[det([[E[x][y] for y in U] for x in S]) for U in subs] for S in subs]
        sP=smith_vals(P,900); sR=smith_vals(P+C,900); sCol=smith_vals([p+c for p,c in zip(P,C)],900)
        print(f"a={a} t={t} z={z} k={k}: row incl C_k in Per_k: {sP==sR}; col incl: {sP==sCol}",flush=True)
