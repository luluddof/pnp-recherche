import numpy as np, sys, warnings
warnings.filterwarnings('ignore')
from trilib import *
fn=sys.argv[1]; tl=float(sys.argv[2]); d=np.load(fn); ps=list(d['ps']); N=4
b=box3([N]*4); pr=Prob(b,ps); x=d['x']
T,R,U,s,V,E,f=tail_state(pr,x)
print(fn,'before tail',f'{f:.3e}','|x|',f'{np.linalg.norm(x):.1f}','sig p0-2..p0+2',np.round(s[ps[0]-3:ps[0]+2],6),flush=True)
x,f,h=lm_solve(pr,x,maxit=100000,tlim=tl)
T,R,U,s,V,E,f=tail_state(pr,x)
print(fn,'after tail',f'{f:.3e}','|x|',f'{np.linalg.norm(x):.1f}','sig',np.round(s[ps[0]-3:ps[0]+2],6),'|A_j|',[float(np.linalg.norm(pr.unpack(x)[2*k]@pr.unpack(x)[2*k+1].T)) for k in range(2)],flush=True)
