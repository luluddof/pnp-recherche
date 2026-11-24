# LP over normalised log-dimension profiles f(E)=L(E)/lambda, E=0..4h:
#  f(0)=0, f nondecreasing, f(h)>=1 (N <= dim V_h), line at every threshold s0: f(2h+s0-1)+f(2h-s0) >= 2 (G20b counting form),
#  cost M >= f(2h+s)+f(2h-s) for |s|<=2h (per-piece, paper's form).  Minimise M (=4 rho).  Optionally concave.
import numpy as np
from scipy.optimize import linprog
def solve(h, concave=False, endpoint_only=False):
    n=4*h+1; nv=n+1  # f_0..f_4h, M
    A=[];b=[]
    def row(): return np.zeros(nv)
    for E in range(n-1):           # monotone f_E - f_{E+1} <= 0
        r=row(); r[E]=1; r[E+1]=-1; A.append(r); b.append(0)
    r=row(); r[h]=-1; A.append(r); b.append(-1)
    rng_s0 = [2*h] if endpoint_only else range(-2*h+1, 2*h+1)
    for s0 in rng_s0:
        r=row(); r[2*h+s0-1]-=1; r[2*h-s0]-=1; A.append(r); b.append(-2)
    for s in range(-2*h,2*h+1):
        r=row(); r[2*h+s]+=1; r[2*h-s]+=1; r[-1]=-1; A.append(r); b.append(0)
    if concave:
        for E in range(1,n-1):
            r=row(); r[E-1]=1; r[E+1]=1; r[E]=-2; A.append(r); b.append(0)
    bounds=[(0,0)]+[(0,None)]*(n-1)+[(None,None)]
    c=np.zeros(nv); c[-1]=1
    res=linprog(c,A_ub=np.array(A),b_ub=np.array(b),bounds=bounds,method='highs')
    return res.fun, res.x[:n]
for h in [5,20,80]:
    M1,f1=solve(h); M2,f2=solve(h,concave=True); M3,_=solve(h,endpoint_only=True); M4,_=solve(h,True,True)
    print('h=%d  no concavity: 4rho>=%.4f | concave: %.4f | endpoint-only floor: %.4f, concave %.4f'%(h,M1,M2,M3,M4))
    if h==20: print(' optimal non-concave profile f at E=0,h/2,h,2h,3h,3.5h,4h:',[round(f1[E],3) for E in [0,10,20,40,60,70,80]])
# step function check: f=0 at 0, 1 on [1,4h-1], 2 at 4h-1..
h=20; f=np.ones(4*h+1); f[0]=0; f[4*h-1:]=2
print('step profile cost max_s:',max(f[2*h+s]+f[2*h-s] for s in range(-2*h,2*h+1)),
      ' line violations:',sum(f[2*h+s0-1]+f[2*h-s0]<2-1e-9 for s0 in range(-2*h+1,2*h+1)))
f=np.ones(4*h+1); f[0]=0
print('step with f(4h)=2 only: cost',max(f[2*h+s]+f[2*h-s] for s in range(-2*h,2*h+1)),' violates endpoint line f(0)+f(4h-1)>=2:',f[0]+f[4*h-1])
