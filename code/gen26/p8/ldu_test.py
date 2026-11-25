# Candidate explicit construction for (A): order rows/cols by (x', rho) (unions of initial segments per class attain Thm-1 minimal minors),
# E = L D U with L unit lower triangular over O (L_ij = ratio of minors, integral since leading minors are minimal).
# F = Per_k(L)^{-1} Per_k(E); row T: w_T = min v2, reduction r_T = row/2^{w_T} mod 2.  Check: w_T <= s(T) for all T and {r_T} independent.
import sys, itertools
from fractions import Fraction
sys.path.insert(0,'.')
from percomp import perm, MOD2
from union_check2 import Eblk
from verify_am2 import v2, v2q
def inv_mod(A):   # inverse of matrix over Z/2^K (must be invertible mod 2)
    n=len(A); M=[r[:]+[int(i==j) for j in range(n)] for i,r in enumerate(A)]
    for c in range(n):
        p=next(r for r in range(c,n) if M[r][c]&1); M[c],M[p]=M[p],M[c]
        iv=pow(M[c][c],-1,MOD2); M[c]=[x*iv%MOD2 for x in M[c]]
        for r in range(n):
            if r!=c and M[r][c]:
                f=M[r][c]; M[r]=[(x-f*y)%MOD2 for x,y in zip(M[r],M[c])]
    return [r[n:] for r in M]
def matmul(A,B): return [[sum(a*b for a,b in zip(r,col))%MOD2 for col in zip(*B)] for r in A]
def run(a,t,zeta,b,kmax):
    E=Eblk(a,t,zeta,b); n=len(E); L_=1<<a; i=v2(t)
    B=[x for x in range(L_) if (x>>i)&1==b]
    order=sorted(range(n),key=lambda k:(B[k]>>(i+1),B[k]%(1<<i)))
    E=[[E[r][c] for c in order] for r in order]
    # LDU without pivoting: L[:,j] = column j of Schur complement / pivot
    A=[r[:] for r in E]; Lm=[[int(r==c) for c in range(n)] for r in range(n)]; s=[]
    for j in range(n):
        piv=A[j][j]; vp=v2(piv); s.append(vp)
        for r in range(j+1,n):
            if A[r][j]:
                assert v2(A[r][j])>=vp, "L not integral"
                f=((A[r][j]>>vp)*pow(piv>>vp,-1,MOD2))%MOD2; Lm[r][j]=f
                A[r]=[(x-f*y)%MOD2 for x,y in zip(A[r],A[j])]
    ok=True
    for k in range(1,min(n,kmax)+1):
        subs=list(itertools.combinations(range(n),k))
        PE=[[perm([[E[x][y] for y in U] for x in S]) for U in subs] for S in subs]
        PL=[[perm([[Lm[x][y] for y in U] for x in S]) for U in subs] for S in subs]
        F=matmul(inv_mod(PL),PE)
        ws=[min(v2(z) for z in row) for row in F]
        sT=[sum(s[j] for j in T) for T in subs]
        le=all(w<=st for w,st in zip(ws,sT))
        basis={}
        for row,w in zip(F,ws):
            v=int(''.join(str((z>>w)&1) for z in row),2)
            while v:
                h=v.bit_length()-1
                if h in basis: v^=basis[h]
                else: basis[h]=v; break
        rk=len(basis)
        print(f"  a={a} t={t} z={zeta} k={k}: w_T<=s(T): {le}  rank {rk}/{len(subs)}  (w_T-s(T) min/max {min(w-st for w,st in zip(ws,sT))}/{max(w-st for w,st in zip(ws,sT))})",flush=True)
        ok&=le and rk==len(subs)
    return ok, s
if __name__=='__main__':
    for (a,t,z) in [(3,1,1),(3,2,1),(3,4,1),(3,-3,Fraction(1,3)),(4,1,1),(4,2,1),(4,4,-3),(4,8,1)]:
        ok,s=run(a,t,z,0,4); print(f"a={a} t={t} z={z}: pivots {s}  construction OK: {ok}",flush=True)
