import time, numpy as np, itertools
d,k=2,4; N=16
T=np.random.rand(N,N,N,N); perms=list(itertools.permutations(range(k)))
t=time.time()
big=T.reshape([d]*(4*k)); acc=np.zeros_like(big)
for pm in perms:
    ax=[]
    for leg in range(4): ax+=[leg*k+pm[i] for i in range(k)]
    acc+=big.transpose(ax)
print("symm",time.time()-t)
t=time.time(); np.linalg.svd(T.reshape(256,256)); print("svd",time.time()-t)
