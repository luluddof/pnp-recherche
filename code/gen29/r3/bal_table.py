# balanced (min max(p,q)) comparison at N=2^k: Majorana formula vs HR tiling vs JH/Hrubes count vs N^2/8
from math import comb, ceil, log2
def rho(N):
    a=0
    while N%2==0: N//=2; a+=1
    c,d=a%4,a//4; return 8*d+2**c
def maj(k):
    n=2*k+1; best=None
    for par in (k%2,):
        um=(n-par)//2
        for lo in range(um+1):
            for hi in range(lo,um+1):
                d=hi-lo
                q=sum(comb(n,j) for j in range(d+1)) if d<=k else 4**k
                p=sum(comb(n,2*u+par) for u in range(um+1) if u<lo or u>hi)
                if best is None or max(p,q)<best: best=max(p,q)
    return best
def jh(r):
    n=N=2**r; t=rho(n)-1; best=N*N
    for k in range(1,t+1):
        if comb(t,k)<N: continue
        for j1 in range(-1,k):
            j2=k-1-j1
            r1=0 if j1<0 else min(comb(t,j1),N); r2=0 if j2<0 else min(comb(t,j2),N)
            best=min(best,max(n*r1,n*r2))
    return best
print(' k   Majorana      JH(Hrubes)    HRtiling      N^2/8   | winner  maj_exp')
for k in range(3,41):
    N=2**k; M=maj(k); J=jh(k); H=N*ceil(N/rho(N)); E=N*N//8
    others=min(J,H,E); w='MAJ' if M<others else ('tie' if M==others else 'other')
    print(f'{k:2d} {M:13d} {J:13d} {H:13d} {E:13d} | {w:5s} {log2(M)/k:.3f} JHexp {log2(J)/k:.3f}')
