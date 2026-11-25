import sys,time; sys.path.insert(0,'.')
from smith7 import *
for (t,zeta,b) in [(1,Fraction(1),0),(-6,Fraction(-3),1),(-(1<<7)*3,Fraction(1,3),0)]:
    t0=time.time(); a=9
    E,i=crossing(a,t,zeta,b); pred=predicted(a,i,zeta); N=max(pred)+20
    got=smith_vals([[red(e,N) for e in row] for row in E],N)
    print("a=9 t",t,"zeta",zeta,"b",b,"i",i,"match",got==pred,"time",round(time.time()-t0,1),flush=True)
