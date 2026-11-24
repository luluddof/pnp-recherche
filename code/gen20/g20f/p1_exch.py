import numpy as np, itertools, math, random
from fractions import Fraction
from scipy.optimize import milp, LinearConstraint, Bounds
def ES(d):  # exact E|S_d|
    return Fraction(sum(math.comb(d,k)*abs(2*k-d) for k in range(d+1)),2**d)
def ESf(d): return Fraction(d*math.comb(d-1,(d-1)//2),2**(d-1))
assert all(ES(d)==ESf(d) for d in range(1,60)); print('E|S_d| closed form exact for d<60')
ok=True; mono=True
for d in range(1,2001):
    e=math.comb(d-1,(d-1)//2)*d/2**(d-1) if d<1000 else math.exp(math.lgamma(d)-math.lgamma((d-1)//2+1)-math.lgamma(d-(d-1)//2)+math.log(d)-(d-1)*math.log(2))
    if not (math.sqrt(d/2)-1e-12<=e<=math.sqrt(d)+1e-12): ok=False
    if d>1 and e/d>prev*(1+1e-12): mono=False
    prev=e/d
print('sqrt(d/2) <= E|S_d| <= sqrt(d) for d<=2000:',ok,'  E|S_d|/d non-increasing:',mono)
def Aint(u,w):
    labs=list(itertools.product((0,1),repeat=u)); idx={l:i for i,l in enumerate(labs)}
    rows=[]
    for a,b in itertools.combinations(labs,2):
        d=sum(x!=y for x,y in zip(a,b))
        for _ in range(w.get(d,0)): rows.append([ (b[i]-a[i]) for i in range(u)])   # delta_e(i)=+-1 if cut (H_i = {lab_i=1}), head b
    D=np.array(rows).T; M=D.shape[1]; best=-1
    for s in itertools.product((-1,1),repeat=u-1):
        s=np.array((1,)+s); A=np.hstack([-(s[:,None]*D),(s[:,None]*D),np.ones((u,1))]); c=np.zeros(2*M+1); c[-1]=-1
        res=milp(c,constraints=[LinearConstraint(A,-np.inf,0),LinearConstraint(np.hstack([np.eye(M),np.eye(M),np.zeros((M,1))]),0,1)],
                 integrality=np.r_[np.ones(2*M),0],bounds=Bounds(np.r_[np.zeros(2*M),-1e4],np.r_[np.ones(2*M),1e4]))
        best=max(best,round(-res.fun))
    pred=sum(Fraction(w.get(d,0)*2**u*math.comb(u,d),2)*ES(d) for d in range(1,u+1))/u
    return best,pred,M
random.seed(2); bad=0
for u in (2,3,4):
    for trial in range(6):
        w={d:random.randint(0,3) for d in range(1,u+1)}
        if sum(w.values())==0: continue
        A,pred,M=Aint(u,w)
        print(f'u={u} w={w} M={M}: A_int={A} formula={pred} ({float(pred):.3f})'); bad+= (A!=pred)
print('mismatches:',bad)
