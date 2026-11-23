import json, itertools, time, numpy as np
from collections import Counter
from mcsp_sat import complexity
# cross-check n=3
b3=json.load(open('c3_brute.json'))
c3=[complexity(f,3) for f in range(256)]
assert c3==b3, "mismatch n=3"
print("n=3 SAT matches brute force:", Counter(c3))
json.dump(c3,open('c3.json','w'))
# n=4 NPN canonical forms
n=4; N=16
F=np.arange(1<<N,dtype=np.int64)
bits=((F[:,None]>>np.arange(N))&1).astype(np.int64)
canon=np.full(1<<N, 1<<20, dtype=np.int64)
w=(1<<np.arange(N)).astype(np.int64)
for perm in itertools.permutations(range(n)):
    for negm in range(1<<n):
        mp=[]
        for t in range(N):
            u=0
            for i in range(n):
                b=((t>>perm[i])&1)^((negm>>i)&1)
                u|=b<<i
            mp.append(u)
        g=bits[:,mp]@w
        canon=np.minimum(canon,np.minimum(g,(1<<N)-1-g))
reps=sorted(set(canon.tolist()))
print("NPN classes:",len(reps))
t0=time.time(); cc={}
for r in reps:
    cc[r]=complexity(r,4)
print("SAT time",time.time()-t0)
c4=[cc[c] for c in canon.tolist()]
print("n=4 distribution over all 65536:",sorted(Counter(c4).items()))
print("class distribution:",sorted(Counter(cc.values()).items()))
json.dump(c4,open('c4.json','w'))
