import numpy as np, pickle, sys, time
from oqgrp import oq_perms
P=oq_perms(3).astype(np.int64)   # (40320,64)
W=np.array([1<<(6*i) for i in range(8)],dtype=np.int64)
def canon(S):
    A=np.sort(P[:,list(S)],axis=1)
    codes=(A*W[:A.shape[1]][::-1]).sum(1)
    i=int(np.argmin(codes)); return tuple(int(t) for t in A[i])
def reps(qq, prev=None):
    R=prev if prev else {(0,)}
    size=len(next(iter(R)))
    while size<qq:
        new=set()
        for S in R:
            for v in range(1,64):
                if v in S: continue
                new.add(canon(S+(v,)))
        R=new; size+=1
        print(size,len(R),flush=True)
    return sorted(R)
if __name__=='__main__':
    q=int(sys.argv[1]); t0=time.time()
    R=reps(q); pickle.dump(R,open(f'reps{q}.pkl','wb')); print(len(R),'reps',time.time()-t0,'s')
