# rank-method (block-compression) lower bound RB(N;p,q) = max over partitions N = sum k_t of sum_t max(0, k N - p - k^2 q, k N - k^2 p - q)
from functools import lru_cache
def RB(N,p,q):
    g=lambda k: max(0, k*N-p-k*k*q, k*N-k*k*p-q)
    @lru_cache(None)
    def best(n): return 0 if n==0 else max(g(k)+best(n-k) for k in range(1,n+1))
    return best(N)
def thmC(N,p,q): return max(0,max(max(k*N-p-k*k*q, k*N-k*k*p-q) for k in range(1,N+1)))
data={(3,1,1):5,(3,2,1):3,(3,2,2):2,(3,3,1):0,(4,1,1):12,(4,2,1):10,(4,3,1):7,(4,4,1):5,(4,5,1):3,(4,2,2):8,(4,3,2):5,(4,5,2):2,(5,1,1):21,(5,2,1):19}
for (N,p,q),d in sorted(data.items()):
    print(N,p,q,'d2=',d,'c=',N*N-d,' thmC(single k)=',thmC(N,p,q),' RB(blocks)=',RB(N,p,q))
