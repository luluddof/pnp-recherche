# Exact check of the Pauli-diagonal dictionary: for f on F_2^{2n}, Phi_f(Z) = sum_a f(a) <W_a,Z> W_a / N.
# Claim: rr(Phi_f) = |supp f| ; rc(Phi_f) = |supp hat f| (omega-Fourier).  tau = Phi_chi.
# Construction tested: AW-type window: f = polynomial P(z) of degree D in z=|a_x & a_z| interpolating (-1)^z on window W;
# g = chi - f (support = classes outside W, plus accidental zeros). R = Phi_g, C = Phi_f.  All exact over Z/Q.
import itertools, sys
from fractions import Fraction
import flint
def run(n,D,c0):
    N=2**n
    pts=list(itertools.product(range(2),repeat=2*n))
    Xm=[[1 if (j^1)==i else 0 for j in range(2)] for i in range(2)]
    def W(a):
        ax,az=a[:n],a[n:]
        M=[[0]*N for _ in range(N)]
        for col in range(N):
            bits=[(col>>(n-1-k))&1 for k in range(n)]
            sign=(-1)**sum(b*z for b,z in zip(bits,az))     # Z^{az}
            row_bits=[b^x for b,x in zip(bits,ax)]            # then X^{ax}
            row=sum(b<<(n-1-k) for k,b in enumerate(row_bits))
            M[row][col]=sign
        return M
    Ws={a:W(a) for a in pts}
    def chi(a): return (-1)**sum(a[k]*a[n+k] for k in range(n))
    # check transpose: W_a^T = chi(a) W_a
    for a in pts:
        M=Ws[a]; assert all(M[j][i]==chi(a)*M[i][j] for i in range(N) for j in range(N))
    z=lambda a: sum(a[k]*a[n+k] for k in range(n))
    win=list(range(c0,c0+D+1))
    # Lagrange interpolation of (-1)^z on win, exact rationals
    def P(t):
        s=Fraction(0)
        for i,xi in enumerate(win):
            term=Fraction((-1)**xi)
            for j,xj in enumerate(win):
                if j!=i: term*=Fraction(t-xj,xi-xj)
            s+=term
        return s
    f={a:P(z(a)) for a in pts}
    g={a:Fraction(chi(a))-f[a] for a in pts}
    om=lambda a,b: sum(a[k]*b[n+k]+a[n+k]*b[k] for k in range(n))%2
    fhat={b:sum(f[a]*(-1)**om(a,b) for a in pts) for b in pts}
    def superop(h):   # matrix S[(x,y),(u,v)] with (Phi_h Z)_{xy} = sum S Z_{uv}
        S=[[Fraction(0)]*(N*N) for _ in range(N*N)]
        for a in pts:
            if h[a]==0: continue
            M=Ws[a]
            for x in range(N):
                for y in range(N):
                    if M[x][y]==0: continue
                    for u in range(N):
                        for v in range(N):
                            if M[u][v]: S[x*N+y][u*N+v]+=h[a]*M[x][y]*M[u][v]/N
        return S
    R=superop(g); C=superop(f)
    # R + C = tau
    for x in range(N):
        for y in range(N):
            for u in range(N):
                for v in range(N):
                    assert R[x*N+y][u*N+v]+C[x*N+y][u*N+v]==(1 if (u==y and v==x) else 0)
    def rank(rows):
        den=1
        for r in rows:
            for e in r: den=den*e.denominator//__import__('math').gcd(den,e.denominator)
        return flint.fmpz_mat([[int(e*den) for e in r] for r in rows]).rank()
    rr=rank(R)
    Cre=[[C[x*N+y][u*N+v] for y in range(N) for v in range(N)] for x in range(N) for u in range(N)]
    rc=rank(Cre)
    sg=sum(1 for a in pts if g[a]!=0); sf=sum(1 for b in pts if fhat[b]!=0)
    print(f"n={n} N={N} D={D} window={win}: rr(R)={rr} |supp g|={sg} ; rc(C)={rc} |supp fhat|={sf} ; R+C=tau OK")
for n,D,c0 in [(1,0,0),(2,0,0),(2,1,0),(2,1,1),(3,0,0),(3,1,0),(3,1,1),(3,2,0),(3,2,1)]:
    run(n,D,c0)
