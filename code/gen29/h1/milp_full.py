# Full Bell-model MILP at small N: choose support (<= q characters, 0 in S WLOG by twisted translation), coefficients
# |kappa_s| <= Cb, maximise exact agreements.  Optimal => exact optimum UNDER the coefficient bound.
import numpy as np, sys, highspy
from bell import *
k=int(sys.argv[1]); q=int(sys.argv[2]); Cb=float(sys.argv[3]); tl=float(sys.argv[4])
n=4**k; b=beta(k); M=charmat(k,list(range(n))).astype(float)
Mbig=Cb*q+1.0
h=highspy.Highs(); h.setOptionValue('output_flag',False); h.setOptionValue('time_limit',tl); h.setOptionValue('threads',1)
inf=highspy.kHighsInf
# vars: kappa (n), used (n), agree (n)
for i in range(n): h.addVar(-Cb,Cb)
for i in range(2*n): h.addVar(0,1)
h.changeColsIntegrality(2*n, np.arange(n,3*n,dtype=np.int32), np.array([highspy.HighsVarType.kInteger]*(2*n)))
h.changeColsCost(3*n, np.arange(3*n,dtype=np.int32), np.concatenate([np.zeros(2*n),-np.ones(n)]))
h.changeColBounds(n, 1.0, 1.0)   # 0 in S
for s in range(n):
    h.addRow(-inf,0,2,np.array([s,n+s],dtype=np.int32),np.array([1.0,-Cb]))
    h.addRow(0,inf,2,np.array([s,n+s],dtype=np.int32),np.array([1.0,Cb]))
h.addRow(-inf,q,n,np.arange(n,2*n,dtype=np.int32),np.ones(n))
for z in range(n):
    idx=np.concatenate([np.arange(n),[2*n+z]]).astype(np.int32)
    h.addRow(-inf, b[z]+Mbig, n+1, idx, np.concatenate([M[z],[Mbig]]))
    h.addRow(b[z]-Mbig, inf, n+1, idx, np.concatenate([M[z],[-Mbig]]))
h.run()
info=h.getInfo(); sol=np.array(h.getSolution().col_value)
x=sol[:n]; x[np.abs(x)<1e-7]=0; agree=int(np.sum(np.abs(M@x-b)<1e-6))
print(f'k={k} q<={q} |kappa|<={Cb}: {h.modelStatusToString(h.getModelStatus())}; found p={n-agree} (nnz {int(np.sum(x!=0))}); '
      f'dual bound p >= {n+info.mip_dual_bound:.2f}; gap {info.mip_gap:.3f}',flush=True)
