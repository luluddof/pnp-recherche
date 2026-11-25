# small-q end of the frontier: delta_N(q) = 2 - log_N f(N,q), delta'_N(q) = log_N( f(N,1)/f(N,q) ) (removes the factor-2 constant
# of f(N,1) = N(N-1)/2), eta = log_N q, for q = 3..N; fits delta' = a*eta (linear corner) vs b*eta^2 (App B shape) vs a*eta+b*eta^2.
import math, numpy as np, json
import analysis as A
rows=[]
for N in range(4,13):
    Fr=A.frontier(N); L=math.log(N); f1=N*(N-1)//2
    def f(q): return min(p for p,qq in Fr if qq<=q)
    for q in range(3,N+1):
        p=f(q)
        if p>=f1: continue
        rows.append((N,q,p,math.log(q)/L,2-math.log(p)/L,math.log(f1/p)/L))
print("N q f(N,q) eta delta delta'")
for r in rows: print(f"{r[0]:2d} {r[1]:2d} {r[2]:3d} {r[3]:.3f} {r[4]:.3f} {r[5]:.3f}")
X=np.array([r[3] for r in rows]); Y=np.array([r[5] for r in rows])
def fit(cols):
    M=np.stack(cols,1); c,res,_,_=np.linalg.lstsq(M,Y,rcond=None); r=Y-M@c; return c, float(np.sqrt(np.mean(r**2)))
for name,cols in [('a*eta',[X]),('b*eta^2',[X**2]),('a*eta+b*eta^2',[X,X**2]),('a*(eta-e0) [free intercept]',[X,np.ones_like(X)])]:
    c,rms=fit(cols); print(f"fit delta' = {name}: coef {np.round(c,3)} rms {rms:.3f}")
# per N: slope of delta' vs eta over q in [N/2, N]
for N in range(6,13):
    sub=[r for r in rows if r[0]==N]
    if len(sub)>=2:
        x=np.array([r[3] for r in sub]); y=np.array([r[5] for r in sub]); c=np.polyfit(x,y,1)
        print(f"N={N}: points {len(sub)}, delta'(eta=1, q=N) = {sub[-1][5]:.3f}, local linear slope {c[0]:.2f}, App-B eta^2/4 at eta=1: 0.25")
json.dump(rows,open('corner_rows.json','w'))
