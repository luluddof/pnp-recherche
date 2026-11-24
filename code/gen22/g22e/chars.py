# exact dimension sequences (integer weights in units 1/den) for several vertex algebras / modules
import math
def prod_series(factors, N):
    # factors: list of (k, e): multiply by (1-q^k)^(-e) if e>0 ; (1+q^k)^m encoded as ('p',k,m)
    a=[0]*(N+1); a[0]=1
    for f in factors:
        if f[0]=='p':
            _,k,m=f
            for _ in range(m):
                for n in range(N,k-1,-1): a[n]+=a[n-k]
        else:
            k,e=f
            for _ in range(e):
                for n in range(k,N+1): a[n]+=a[n-k]
    return a
def mul(a,b,N):
    c=[0]*(N+1)
    for i,x in enumerate(a[:N+1]):
        if x:
            for j in range(0,N+1-i): c[i+j]+=x*b[j]
    return c
def sigma3(n): return sum(d**3 for d in range(1,n+1) if n%d==0)
def E4(N): return [1]+[240*sigma3(n) for n in range(1,N+1)]
def heis(N,d=1): return prod_series([(k,d) for k in range(1,N+1)],N)
def eta_inv_pow(N,d): return heis(N,d)
def J_dims(N):
    # V natural: dims V_n = c(n-1) of J = E4^3/Delta - 744 ; Delta = q prod(1-q^n)^24
    M=N+2
    e=E4(M); e3=mul(mul(e,e,M),e,M)
    inv=eta_inv_pow(M,24)        # prod (1-q^n)^-24
    j=mul(e3,inv,M)              # q*J(q) +? : E4^3/Delta = q^-1 * j  -> j[n] = coeff q^{n-1} of E4^3/Delta
    d=[j[n] for n in range(0,N+1)]  # coefficient of q^{n-1} in J+744 -> dim V_n for n!=1
    d[1]-=744
    return d
def leech(N):
    d=J_dims(N); d[1]+=24; return d
def E8(N): return mul(E4(N),heis(N,8),N)
def NSferm(N2,n=1):  # grading in half units, length N2+1
    return prod_series([('p',2*k-1,n) for k in range(1,N2//2+2)],N2)
def ramond(N,n=1):   # integer grading; factor 2^{n/2} handled separately
    return prod_series([('p',k,n) for k in range(1,N+1)],N)
def RR(N,res,mod=5):
    return prod_series([(k,1) for k in range(1,N+1) if k%mod in res],N)
def symplf(N,n=1): return prod_series([('p',k,2*n) for k in range(1,N+1)],N)
