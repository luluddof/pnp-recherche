import sys, math, numpy as np, importlib.util
spec=importlib.util.spec_from_file_location('appb','/home/user/pnp-recherche/code/gen20/g20b/appb_exact.py'); appb=importlib.util.module_from_spec(spec); spec.loader.exec_module(appb)
import gauge
h=int(sys.argv[1]); N,mus,I,T=appb.pieces(h); d=math.factorial(h); n=N*N
mat=lambda T4: np.array([[T4[x][y][u][v] for u in range(N) for v in range(N)] for x in range(N) for y in range(N)],dtype=float)
S=np.zeros((n,n))
for x in range(N):
    for y in range(N): S[x*N+y,y*N+x]=1
R=np.zeros((n,n))
for s0 in I:
    C=d*S-R; K=C/d
    Kr=K.reshape(N,N,N,N).transpose(0,2,1,3).reshape(n,n); U,sg,Vt=np.linalg.svd(Kr); q=int(np.sum(sg>1e-9*sg[0]))
    p=np.linalg.matrix_rank(S-K,tol=1e-9)
    if 0<p<n and 0<q<n and abs(p-q)<=max(p,q):
        G=np.array([np.sqrt(sg[i])*U[:,i].reshape(N,N) for i in range(q)]); H=np.array([np.sqrt(sg[i])*Vt[i].reshape(N,N) for i in range(q)])
        assert np.allclose(np.einsum('kxu,kyv->xyuv',G,H).reshape(n,n),K)
        Q0=gauge.Qop(G,H); G2,H2=gauge.run(G,H,N,1); Q1=gauge.Qop(G2,H2); t1,_=gauge.Tsnorm(G2,H2,N)
        print(f"h={h} N={N} s0={s0} p={p} q={q}  Q(svd)={Q0:.3g}  Q after X,Y gauge={Q1:.3g}  ||Ts|| after={t1:.3g}",flush=True)
    R=R+mat(T[s0])
