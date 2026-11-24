# H = 0 example: commuting tensor factors, labels in top level W1(x)W2 (dims 2,2), o(a_x)=A_x (x)1, o(b_v)=1(x)B_v.
# find A_x, B_v (2x2), labels Y_y, U_u in C^2(x)C^2 with <U_u|A_x(x)B_v|Y_y> = [x=v][y=u], <U_u|Y_y> = [u=y], N = 2.
import numpy as np, scipy.optimize as so
rng=np.random.default_rng(1); N=2
def unpack(z):
    A=z[:8].reshape(2,2,2); B=z[8:16].reshape(2,2,2); Y=z[16:24].reshape(2,4); U=z[24:32].reshape(2,4); return A,B,Y,U
def res(z):
    A,B,Y,U=unpack(z); r=[]
    for x in range(N):
        for v in range(N):
            K=np.kron(A[x],B[v])
            for u in range(N):
                for y in range(N): r.append(U[u]@K@Y[y]-(x==v)*(u==y))
    for u in range(N):
        for y in range(N): r.append(U[u]@Y[y]-(u==y))
    return np.array(r)
best=None
for t in range(200):
    s=so.least_squares(res,rng.normal(size=32),xtol=1e-15,ftol=1e-15,gtol=1e-15)
    if best is None or s.cost<best.cost: best=s
    if s.cost<1e-24: break
print("restarts",t+1,"residual",np.abs(res(best.x)).max())
A,B,Y,U=unpack(best.x); np.set_printoptions(precision=4,suppress=True)
print("A",A,"\nB",B,"\nY",Y,"\nU",U)
