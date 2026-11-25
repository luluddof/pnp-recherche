# Apply verified certificates to the 5 exact optima: check tr g(Ts/Q) in [a m-(n-m), L N] and bound <= p.
import json, glob, numpy as np
from fractions import Fraction as F
from common import *
certs=[json.load(open(f)) for f in sorted(glob.glob('cert_*.json'))]
for f,pst in OPT:
    N,A,Bm=load(f); q=len(A); n=N*N; Th=Kmat(A,Bm)@swap(N); Ts=(Th+Th.T)/2
    p=np.linalg.matrix_rank(np.eye(n)-Th); m=n-p; ev=np.linalg.eigvalsh(Ts)
    for c in certs:
        Q=c['q']
        if Q<q: continue
        co=[float(F(x)) for x in c['coeffs']]; a=float(F(c['a'])); L=float(F(c['L']))
        tr=sum(np.polyval(co[::-1],ev/Q))
        lo=a*m-(n-m); hi=L*N; bound=(a*n-L*N)/(1+a)
        print(f'{f.split("/")[-1]} q={q} cert(Q={Q},a={a}): lower {lo:.2f} <= tr g = {tr:.2f} <= {hi:.2f} : {lo-1e-9<=tr<=hi+1e-9};  bound {bound:.2f} <= p={p}: {bound<=p}')
