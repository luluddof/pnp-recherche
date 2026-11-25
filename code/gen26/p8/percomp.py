# Test (A): Smith of the permanental compound Per_k(E) (rows S, cols U, |S|=|U|=k, entry per E[S,U]) vs k-subset sums of Smith(E).
# (A) is necessary for the graded strengthening C1'' of C1'-right. Also report full N=M(E) block structure check.
import sys, itertools
sys.path.insert(0,'.')
sys.path.insert(0,'/home/user/pnp-recherche/code/gen26/p6')
from verify_am2 import v2, MOD
from core_explore import Ecore
exec(open('/home/user/pnp-recherche/code/gen26/p6/m_profile.py').read().split('def Eblock')[0].split('K = 700')[1].replace('; MOD = 1 << K','',1) if False else '')
K2=700; MOD2=1<<K2
def smith_vals(A, cap):
    A=[r[:] for r in A]; rows=set(range(len(A))); cols=list(range(len(A[0]))); out=[]
    while rows and cols:
        best=None
        for r in rows:
            Ar=A[r]
            for c in cols:
                e=Ar[c]
                if e:
                    v=(e&-e).bit_length()-1
                    if best is None or v<best[0]: best=(v,r,c)
            if best and out and best[0]==out[-1]: break
        if best is None or best[0]>=cap: break
        v,pr,pc=best; out.append(v); inv=pow(A[pr][pc]>>v,-1,MOD2); P=A[pr]
        nz=[c for c in cols if P[c]]
        for r in rows:
            if r!=pr and A[r][pc]:
                m=((A[r][pc]>>v)*inv)%MOD2; Ar=A[r]
                for c in nz: Ar[c]=(Ar[c]-m*P[c])%MOD2
        rows.discard(pr); cols.remove(pc)
    return sorted(out)
def perm(M):
    n=len(M)
    if n==0: return 1
    # Ryser
    tot=0
    for mask in range(1,1<<n):
        cols=[j for j in range(n) if mask>>j&1]
        p=1
        for r in M: p=p*sum(r[j] for j in cols)%MOD2
        tot+=(-1)**len(cols)*p
    return (tot*(-1)**n)%MOD2
def sform(a,i,v2z=0): return [(a-i+1)+v2z+2*k+2*sum(v2(m) for m in range(1,k+1)) for k in range(1<<(a-i-1))]*(1<<i)
if __name__=='__main__':
    a=int(sys.argv[1]); i=int(sys.argv[2]); kmax=int(sys.argv[3])
    E,_=Ecore(a,1<<i,1); n=len(E); E=[[x%MOD2 for x in r] for r in E]
    s=sform(a,i)
    for k in range(1,kmax+1):
        subs=list(itertools.combinations(range(n),k))
        P=[[perm([[E[x][y] for y in U] for x in S]) for U in subs] for S in subs]
        sm=smith_vals(P,600)
        sq=sorted(sum(c) for c in itertools.combinations(s,k))
        ok=all(p<=q for p,q in zip(sm,sq))
        print(f"a={a} i={i} k={k}: Smith(Per_k)<=k-sqfree termwise: {ok}\n   Per: {sm[:40]}\n   sqf: {sq[:40]}",flush=True)
