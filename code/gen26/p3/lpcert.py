# Best bound obtainable by the odd-trace method (Theorems A/B):  p >= (a n - N L)/(1+a),
# L = min sum_j |g_j| over odd polys g = sum g_j y^j (deg <= k) with g <= 1 on [0,1], g >= a on [1/q,1].
import numpy as np, sys
from scipy.optimize import linprog
def Lstar(q,a,k,grid=1500):
    js=list(range(1,k+1,2)); J=len(js)
    ys=np.concatenate([np.linspace(0,1,grid),np.geomspace(1e-4,1,grid)]); ys=np.unique(ys)
    V=np.array([[y**j for j in js] for y in ys])
    # vars: g (J), u (J) ; min sum u ; -u<=g<=u ; V g <= 1 ; -V g <= -a on y>=1/q
    c=np.concatenate([np.zeros(J),np.ones(J)])
    A=[];b=[]
    I=np.eye(J)
    A.append(np.hstack([I,-I])); b.append(np.zeros(J))
    A.append(np.hstack([-I,-I])); b.append(np.zeros(J))
    A.append(np.hstack([V,np.zeros((len(ys),J))])); b.append(np.ones(len(ys)))
    sel=ys>=1.0/q
    A.append(np.hstack([-V[sel],np.zeros((sel.sum(),J))])); b.append(-a*np.ones(sel.sum()))
    r=linprog(c,A_ub=np.vstack(A),b_ub=np.concatenate(b),bounds=[(None,None)]*J+[(0,None)]*J,method='highs')
    return r.fun if r.status==0 else np.inf
def best(N,q,kmax=41):
    n=N*N; out=(0,None)
    for k in range(1,kmax+1,4):
        for a in [1/q*t for t in (1,1.5,2,3,4,6,8,12,16,24,32)]+[0.25,0.5,0.75,0.9,1.0]:
            L=Lstar(q,a,k)
            if L<np.inf:
                val=(a*n-N*L)/(1+a)
                if val>out[0]: out=(val,(k,round(a,4),round(L,2)))
    return out
if __name__=='__main__':
    for q in (3,10,30):
        for N in (10**3,10**6,10**9,10**12):
            v,arg=best(N,q)
            print('q',q,'N %.0e'%N,'best odd-trace bound p >= %.4g N^2/q'%(v*q/N**2),' (p/N^2=%.4f)'%(v/N**2),'(k,a,L)=',arg,'  ln N/q=%.2f'%(np.log(N)/q),flush=True)
