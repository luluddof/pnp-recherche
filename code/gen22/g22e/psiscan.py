import math, numpy as np
from chars import *
from scaling import mk_phi, run
def scan(name,d,den=1,mult0=1.0):
    phi,mean=mk_phi(d,den)
    out=[]
    for b in [0.2,0.5,1,2,3,4,6,8,12]:
        m,tail=mean(b)
        if tail>1e-12: out.append((b,None)); continue
        out.append((b,round(phi(b)+math.log(mult0)-b*m,4)))
    print("%-28s psi(b)=phi+b phi':"%name, out)
N=400
nsf=NSferm(2*N)
ising_vac=[nsf[2*k] for k in range(N+1)]          # integer-weight part
ising_eps=[nsf[2*k+1] for k in range(N)]          # half-integer part, from its bottom
for nm,d,den,m in [("Ising vacuum",ising_vac,1,1),("Ising eps",ising_eps,1,1),("Ising sigma",ramond(N),1,1),
                   ("Ramond 1 fermion (x sqrt2)",ramond(N),1,math.sqrt(2)),("NS fermion",nsf,2,1),("Heisenberg",heis(N),1,1),
                   ("Yang-Lee vac",RR(N,{2,3}),1,1),("Yang-Lee h=-1/5",RR(N,{1,4}),1,1),("Moonshine",J_dims(300),1,1),("E8",E8(N),1,1)]:
    scan(nm,d,den,m)
for nm,d in [("Ising vacuum",ising_vac),("Ising eps",ising_eps),("Ising sigma",ramond(N))]:
    run(nm+" (scaling)",d)
