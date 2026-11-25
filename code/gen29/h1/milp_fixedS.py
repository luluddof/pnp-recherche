# MILP: max #exact agreements for FIXED support S (Majorana degree<=d frame), |kappa_s| <= Cb, big-M formulation (HiGHS).
import numpy as np, sys, itertools, highspy
from bell import *
k=int(sys.argv[1]); d=int(sys.argv[2]); Cb=float(sys.argv[3]); tl=float(sys.argv[4])
n=4**k; b=beta(k); msk=(1<<k)-1
A=[]
for j in range(k):
    for P in 'XY':
        a=1<<j; bb=((1<<j)-1) | ((1<<j) if P=='Y' else 0); A.append(a|(bb<<k))
t=0
for v in A: t^=v
A.append(t)
def Jm(v): return ((v&msk)<<k) | (v>>k)
S=set([0])
for r in range(1,d+1):
    for T in itertools.combinations(range(2*k+1),r):
        v=0
        for i in T: v^=Jm(A[i])
        S.add(v)
S=sorted(S); q=len(S); M=charmat(k,S).astype(float)
Mbig=Cb*q+1.0
h=highspy.Highs(); h.setOptionValue('output_flag',False); h.setOptionValue('time_limit',tl)
inf=highspy.kHighsInf
for i in range(q): h.addVar(-Cb,Cb)
for z in range(n): h.addVar(0,1)
h.changeColsIntegrality(n, np.arange(q,q+n,dtype=np.int32), np.array([highspy.HighsVarType.kInteger]*n))
h.changeColsCost(q+n, np.arange(q+n,dtype=np.int32), np.concatenate([np.zeros(q),-np.ones(n)]))
for z in range(n):
    idx=np.concatenate([np.arange(q),[q+z]]).astype(np.int32)
    # M x - b <= Mbig(1-a)  ->  M x + Mbig a <= b + Mbig
    h.addRow(-inf, b[z]+Mbig, q+1, idx, np.concatenate([M[z],[Mbig]]))
    # M x - b >= -Mbig(1-a) -> M x - Mbig a >= b - Mbig
    h.addRow(b[z]-Mbig, inf, q+1, idx, np.concatenate([M[z],[-Mbig]]))
h.run()
info=h.getInfo(); sol=h.getSolution().col_value
x=np.array(sol[:q]); agree=int(np.sum(np.abs(M@x-b)<1e-6))
print(f'k={k} d={d} |S|={q} bound |kappa|<={Cb}: status {h.modelStatusToString(h.getModelStatus())}, best agree {agree} (p={n-agree}), '
      f'dual bound p >= {n+info.mip_dual_bound:.1f}, gap {info.mip_gap:.3f}')
