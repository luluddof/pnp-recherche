# Exact p_min(q) at N=8 via: O^+(6,2)-orbit reps of supports (0 in S), and per rep a fixed-support MILP with the RIGOROUS
# coefficient bound |kappa| <= Hadamard(q)/2^{q-1} (Cramer on a q x q +-1 system; det of +-1 q x q matrix is a multiple of 2^{q-1}).
import numpy as np, sys, time, pickle, os, highspy
from bell import beta, charmat
import exact8
q=int(sys.argv[1]); part=int(sys.argv[2]); nparts=int(sys.argv[3]); tl=float(sys.argv[4])
k=3; n=64; b=beta(k)
fn=f'reps{q}.pkl'
cut=float(sys.argv[5]) if len(sys.argv)>5 else None
if os.path.exists(fn): RS=pickle.load(open(fn,'rb'))
else:
    t0=time.time(); RS=exact8.reps(q); pickle.dump(RS,open(fn,'wb')); print(f'{len(RS)} reps in {time.time()-t0:.0f}s',flush=True)
Cb=q**(q/2)/2**(q-1)
def solve(S):
    M=charmat(k,list(S)).astype(float); Mbig=Cb*q+1.0
    h=highspy.Highs(); h.setOptionValue('output_flag',False); h.setOptionValue('time_limit',tl); h.setOptionValue('threads',1)
    if cut is not None: h.setOptionValue('objective_bound', -(n-cut)+0.5)
    inf=highspy.kHighsInf
    for i in range(q): h.addVar(-Cb,Cb)
    for z in range(n): h.addVar(0,1)
    h.changeColsIntegrality(n, np.arange(q,q+n,dtype=np.int32), np.array([highspy.HighsVarType.kInteger]*n))
    h.changeColsCost(q+n, np.arange(q+n,dtype=np.int32), np.concatenate([np.zeros(q),-np.ones(n)]))
    for z in range(n):
        idx=np.concatenate([np.arange(q),[q+z]]).astype(np.int32)
        h.addRow(-inf, b[z]+Mbig, q+1, idx, np.concatenate([M[z],[Mbig]]))
        h.addRow(b[z]-Mbig, inf, q+1, idx, np.concatenate([M[z],[-Mbig]]))
    h.run(); info=h.getInfo(); x=np.array(h.getSolution().col_value[:q])
    agree=int(np.sum(np.abs(M@x-b)<1e-6)) if len(x)==q else 0
    return n-agree, n+info.mip_dual_bound, h.modelStatusToString(h.getModelStatus())
best=None; worst_bound=None; t0=time.time(); notopt=0
for S in RS[part::nparts]:
    p,lb,st=solve(S)
    if st not in ('Optimal','Infeasible'): notopt+=1; print('  nonopt',st,[int(t) for t in S],flush=True)
    if best is None or p<best[0]: best=(p,[int(t) for t in S])
    lbc=np.ceil(lb-1e-6) if st=='Optimal' else (cut+1 if st=='Infeasible' else -1)
    if worst_bound is None or lbc<worst_bound: worst_bound=lbc
print(f'q={q} part {part}/{nparts}: {len(RS[part::nparts])} reps, best p={best}, min over reps of MILP lower bound = {worst_bound}, '
      f'non-optimal solves {notopt}, {time.time()-t0:.0f}s',flush=True)
