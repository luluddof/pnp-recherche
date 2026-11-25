# gauge-minimise Q = sum_l ||X G_l Y^-1|| ||Y H_l X^-1|| (op norms) over X,Y in GL_N, via Frobenius proxy then report op-norm Q,
# ||Ts|| after gauge, and the Theta spectrum (gauge invariant).  usage: gauge.py z.npy N q
import numpy as np, sys
from scipy.optimize import minimize
from scipy.linalg import expm
def load(f,N,q):
    z=np.load(f); s=N*N; return z[:q*s].reshape(q,N,N), z[q*s:].reshape(q,N,N)
def apply(G,H,a,b,N):
    X=expm(a.reshape(N,N)); Y=expm(b.reshape(N,N)); Xi=np.linalg.inv(X); Yi=np.linalg.inv(Y)
    return np.array([X@g@Yi for g in G]), np.array([Y@h@Xi for h in H])
def Qop(G,H): return sum(np.linalg.norm(g,2)*np.linalg.norm(h,2) for g,h in zip(G,H))
def Qfro(G,H): return sum(np.linalg.norm(g)*np.linalg.norm(h) for g,h in zip(G,H))
def Tsnorm(G,H,N):
    s=N*N; K=np.einsum('kxu,kyv->xyuv',G,H).reshape(s,s); S=np.zeros((s,s))
    for x in range(N):
        for y in range(N): S[x*N+y,y*N+x]=1
    Th=K@S; Ts=(Th+Th.T)/2; ev=np.linalg.eigvalsh(Ts); return abs(ev).max(), Th
def run(G,H,N,restarts=2):
    I=np.eye(N)
    def f(w):
        X=I+w[:N*N].reshape(N,N); Y=I+w[N*N:].reshape(N,N)
        try: Xi=np.linalg.inv(X); Yi=np.linalg.inv(Y)
        except np.linalg.LinAlgError: return 1e300, np.zeros_like(w)
        gX=np.zeros((N,N)); gY=np.zeros((N,N)); val=0.
        for g,h in zip(G,H):
            A=X@g@Yi; B=Y@h@Xi; val+=np.sum(A*A)+np.sum(B*B)
            gX+=2*A@(g@Yi).T - 2*Xi.T@B.T@B@Xi.T*0  # placeholder replaced below
            gX-= 0
        return val, None
    def fg(w):
        X=I+w[:N*N].reshape(N,N); Y=I+w[N*N:].reshape(N,N)
        if abs(np.linalg.det(X))<1e-300 or abs(np.linalg.det(Y))<1e-300: return 1e300, np.zeros_like(w)
        Xi=np.linalg.inv(X); Yi=np.linalg.inv(Y)
        gX=np.zeros((N,N)); gY=np.zeros((N,N)); val=0.
        for g,h in zip(G,H):
            A=X@g@Yi; B=Y@h@Xi; val+=np.sum(A*A)+np.sum(B*B)
            # d||X g Yi||^2: dX -> 2 A (g Yi)^T ; dY -> -2 Yi^T A^T A Yi^T ... computed as -(Yi)^T (A^T A) Yi^T? use M=A: d(Yi)=-Yi dY Yi
            gX+=2*A@(g@Yi).T; gY+=-2*(X@g@Yi).T@A@Yi.T*0 + (-2*(Yi.T@ (X@g).T @ A @ Yi.T)).T*0
            gY+=-2*Yi.T@(A.T@ (X@g)).T*0
            # clean formulas: d||A||^2 wrt Y = -2 (Yi^T) (A^T X g)^T ... do it numerically-safe via trace identities:
        return val, None
    def fg2(w):
        X=I+w[:N*N].reshape(N,N); Y=I+w[N*N:].reshape(N,N)
        Xi=np.linalg.inv(X); Yi=np.linalg.inv(Y)
        gX=np.zeros((N,N)); gY=np.zeros((N,N)); val=0.
        for g,h in zip(G,H):
            A=X@g@Yi; B=Y@h@Xi; val+=np.sum(A*A)+np.sum(B*B)
            # ||A||^2, A = X g Yi: dA = dX g Yi - X g Yi dY Yi = dX g Yi - A dY Yi
            gX+=2*A@(g@Yi).T;  gY+=-2*A.T@A@Yi.T
            # ||B||^2, B = Y h Xi: dB = dY h Xi - B dX Xi
            gY+=2*B@(h@Xi).T;  gX+=-2*B.T@B@Xi.T
        return np.log(val), np.concatenate([gX.ravel(),gY.ravel()])/val
    best=None
    rng=np.random.default_rng(0)
    for r in range(restarts):
        w0=np.zeros(2*N*N) if r==0 else rng.standard_normal(2*N*N)*0.05
        res=minimize(fg2,w0,jac=True,method='L-BFGS-B',options=dict(maxiter=5000))
        X=I+res.x[:N*N].reshape(N,N); Y=I+res.x[N*N:].reshape(N,N); Xi=np.linalg.inv(X); Yi=np.linalg.inv(Y)
        G2=np.array([X@g@Yi for g in G]); H2=np.array([Y@h@Xi for h in H]); Q=Qop(G2,H2)
        if best is None or Q<best[0]: best=(Q,G2,H2)
    return best[1],best[2]
if __name__=='__main__':
    f,N,q=sys.argv[1],int(sys.argv[2]),int(sys.argv[3]); G,H=load(f,N,q)
    # per-term balance first
    sc=[np.sqrt(np.linalg.norm(h)/np.linalg.norm(g)) for g,h in zip(G,H)]; G=np.array([g*c for g,c in zip(G,sc)]); H=np.array([h/c for h,c in zip(H,sc)])
    t0,Th=Tsnorm(G,H,N); print(f"{f}: before: Q={Qop(G,H):.3g} ||Ts||={t0:.3g} maxcoef={max(abs(G).max(),abs(H).max()):.3g}")
    G2,H2=run(G,H,N); t1,Th2=Tsnorm(G2,H2,N)
    ev=np.linalg.eigvals(Th2); evs=sorted(ev,key=lambda c:-abs(c))
    print(f"   after gauge: Q={Qop(G2,H2):.3g} ||Ts||={t1:.3g} maxcoef={max(abs(G2).max(),abs(H2).max()):.3g}; Theta top eig {np.round(evs[:5],2)}; #|eig-1|<1e-2: {int(np.sum(abs(ev-1)<1e-2))}")
    np.save(f.replace('.npy','_g.npy'),np.concatenate([G2.ravel(),H2.ravel()]))
