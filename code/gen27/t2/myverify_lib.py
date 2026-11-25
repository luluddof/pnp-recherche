# Independent exact verifier (G26-P15). Pure Python integers / Fractions, no numpy/flint.
import json, sys, glob
from fractions import Fraction

def rank_exact(M):
    # Fraction Gaussian elimination
    A=[[Fraction(v) for v in r] for r in M]; m=len(A); n=len(A[0]); r=0
    for c in range(n):
        piv=next((i for i in range(r,m) if A[i][c]!=0),None)
        if piv is None: continue
        A[r],A[piv]=A[piv],A[r]
        for i in range(r+1,m):
            if A[i][c]!=0:
                f=A[i][c]/A[r][c]
                A[i]=[a-f*b for a,b in zip(A[i],A[r])]
        r+=1
        if r==m: break
    return r

def rank_bareiss(M):
    A=[[int(v) for v in r] for r in M]; m=len(A); n=len(A[0]); r=0; prev=1
    for c in range(n):
        piv=next((i for i in range(r,m) if A[i][c]!=0),None)
        if piv is None: continue
        A[r],A[piv]=A[piv],A[r]
        for i in range(r+1,m):
            A[i]=[(A[r][c]*A[i][j]-A[i][c]*A[r][j])//prev for j in range(n)]
        prev=A[r][c]; r+=1
        if r==m: break
    return r

def check(fn,p,q):
    d=json.load(open(fn)); N=d["N"]; E=d["E"]; H=d["C"]
    for T in E+H:
        assert len(T)==N and all(len(row)==N for row in T)
        for row in T:
            for v in row: assert isinstance(v,int), v
    nterms=len(E); assert len(H)==nterms
    s=N*N
    # rows (x,y), cols (u,v): X = [x=v][y=u]; C[(x,y),(u,v)] = sum_l E_l[x,u] H_l[y,v]
    M=[[0]*s for _ in range(s)]
    for x in range(N):
        for y in range(N):
            row=M[x*N+y]
            for u in range(N):
                for v in range(N):
                    val=(1 if (x==v and y==u) else 0)
                    for l in range(nterms):
                        a=E[l][x][u]
                        if a: val-=a*H[l][y][v]
                    row[u*N+v]=val
    r1=rank_exact(M); r2=rank_bareiss(M)
    # operator-Schmidt rank of C: matrix rows (x,u), cols (y,v)
    R=[[sum(E[l][x][u]*H[l][y][v] for l in range(nterms)) for y in range(N) for v in range(N)] for x in range(N) for u in range(N)]
    osr=rank_bareiss(R)
    # sanity: swap itself has rank N^2, and C=0 must give N^2
    return N,r1,r2,nterms,osr

