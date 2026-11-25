# Sign ansatz for the Z2^k orbit family: X in {+-1}^{NxN}, Y = X^T (scaled 1/2): block kappa = I - Z_kappa,
# Z_kappa[i,j] = X[i,j^k] X[j,i^k]/2 ; want Z_kappa = I - v v^T/2  <=>  X[i,j^k]X[j,i^k] = -X[0,i^k]X[i,k]X[0,j^k]X[j,k] (0<i<j)
# (and diag automatically +1/2).  Linear over F2 in x = (1-X)/2.  Then p = N, q <= N.
import numpy as np, sys, itertools
def solve_f2(rows,rhs,nv):
    # Gaussian elimination over F2 with python ints as bitsets
    piv={}  # pivot col -> (row bits, rhs)
    for r,b in zip(rows,rhs):
        while r:
            c=r.bit_length()-1
            if c in piv: pr,pb=piv[c]; r^=pr; b^=pb
            else: piv[c]=(r,b); break
        else:
            if b: return None,len(piv)
    # back-substitute: free vars = 0
    x=[0]*nv
    for c in sorted(piv):
        r,b=piv[c]; v=b; rr=r^(1<<c)
        while rr:
            cc=rr.bit_length()-1; v^=x[cc]; rr^=(1<<cc)
        x[c]=v
    return x,len(piv)
def system(N):
    var=lambda i,u: i*N+u
    rows=[];rhs=[]
    for k in range(N):
        for i in range(1,N):
            for j in range(i+1,N):
                vs=[var(i,j^k),var(j,i^k),var(0,i^k),var(i,k),var(0,j^k),var(j,k)]
                r=0
                for v in vs: r^=(1<<v)
                rows.append(r); rhs.append(1)
    return rows,rhs
if __name__=='__main__':
    for N in map(int,sys.argv[1:]):
        rows,rhs=system(N); x,rk=solve_f2(rows,rhs,N*N)
        print(f"N={N}: eqs={len(rows)} vars={N*N} rank={rk}:", "SOLVABLE" if x is not None else "INCONSISTENT",flush=True)
        if x is not None:
            X=1-2*np.array(x).reshape(N,N); np.save(f"signX_{N}.npy",X)
