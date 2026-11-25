import sys,random,time
sys.path.insert(0,'.')
from smith7 import *
random.seed(11)
zetas=[Fraction(1),Fraction(-3),Fraction(1,3),Fraction(2),Fraction(-12),Fraction(5,7)]
for a in (6,7):
    L=1<<a; nf=0; nc=0; t0=time.time()
    ts=[t for t in range(-L,L) if t!=0 and v2(t)<a]
    ts=sorted(set(random.sample(ts,10)+[-L+1,L-1,-(L//2),3*(L//4),-(L//4)*3,1,-1,2,-2,5*(L//8)]))
    for t in ts:
        for zeta in random.sample(zetas,2):
            for b in (0,1):
                E,i=crossing(a,t,zeta,b); pred=predicted(a,i,zeta); N=max(pred)+20
                got=smith_vals([[red(e,N) for e in row] for row in E],N); nc+=1
                if got!=pred: nf+=1; print("FAIL",a,t,zeta,b,got,pred,flush=True)
    print("a",a,"cases",nc,"fails",nf,"time",time.time()-t0,flush=True)
