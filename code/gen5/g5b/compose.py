# End-to-end: build Hurwitz-type identity h!(sum a^2)(sum b^2) = sum_k f_k g_k from pieces, verify symbolically.
import sys, io, contextlib, sympy as sp
import indep
h=int(sys.argv[1])
with contextlib.redirect_stdout(io.StringIO()): A,P=indep.run(h)
M=len(P); a=sp.symbols(f'a0:{M}'); b=sp.symbols(f'b0:{M}')
terms=0; total=0
for s,T in A.items():
    R=sp.Matrix([[T[x][y][u][v] for u in range(M) for v in range(M)] for x in range(M) for y in range(M)])
    C=sp.Matrix([[T[x][y][u][v] for y in range(M) for v in range(M)] for x in range(M) for u in range(M)])
    if R.rank()<=C.rank():
        L,Rt=R.rank_decomposition()   # R = L*Rt ; L rows (x,y), Rt cols (u,v)
        for k in range(L.shape[1]):
            f=sum(L[x*M+y,k]*a[x]*b[y] for x in range(M) for y in range(M))
            g=sum(Rt[k,u*M+v]*a[v]*b[u] for u in range(M) for v in range(M))
            total+=f*g; terms+=1
    else:
        L,Rt=C.rank_decomposition()   # rows (x,u), cols (y,v)
        for k in range(L.shape[1]):
            f=sum(L[x*M+u,k]*a[x]*b[u] for x in range(M) for u in range(M))
            g=sum(Rt[k,y*M+v]*a[v]*b[y] for y in range(M) for v in range(M))
            total+=f*g; terms+=1
lhs=sp.factorial(h)*sum(t**2 for t in a)*sum(t**2 for t in b)
print(f"h={h} N={M}: products used={terms}; identity holds: {sp.expand(total-lhs)==0}")
