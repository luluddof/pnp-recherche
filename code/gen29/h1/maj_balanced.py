# Balanced (max(p,q) minimal) points of the Majorana family by exact formula (q = sum_{j<=d} C(2k+1,j), upper bound,
# equal to the WHT count for k<=7), vs singlet family and Hurwitz-Radon tiling s = N*ceil(N/rho(N)).
from math import comb, log2, ceil
def rho(N):
    a=0
    while N%2==0: N//=2; a+=1
    c,dd=divmod(a,4); return 2**dd+8*c
def maj(k):
    n=2*k+1; par=k%2  # parity of J: k=3 ->1, k=4->0 (empirical: par = k mod 2)
    best=None
    for lo in range(k+1):
        for hi in range(lo,k+1):
            d=hi-lo; q=sum(comb(n,j) for j in range(min(d,k)+1)) if d<=k else 4**k
            p=sum(comb(n,2*u+par) for u in range(k+1) if (u<lo or u>hi) and 2*u+par<=n)
            m=max(p,q)
            if best is None or m<best[0]: best=(m,p,q,lo,hi)
    return best
def sing(k):
    best=None
    for lo in range(k+1):
        for hi in range(lo,k+1):
            d=hi-lo; q=sum(comb(k,j)*3**j for j in range(d+1))
            p=sum(comb(k,m)*3**(k-m) for m in range(k+1) if m<lo or m>hi)
            mm=max(p,q)
            if best is None or mm<best[0]: best=(mm,p,q)
    return best
for k in list(range(3,16))+[20,24,31]:
    N=2**k; M=maj(k); S=sing(k); hr=N*ceil(N/rho(N))
    print(f'k={k:2d} N={N}: Majorana bal {M[0]} (p={M[1]},q={M[2]}) exp {log2(M[0])/k:.3f} | singlet {S[0]} exp {log2(S[0])/k:.3f} | HR {hr} | N^2/8 {N*N//8}')
