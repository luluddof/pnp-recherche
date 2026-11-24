"""G8c check: (a) eight-square identity (12) with sigma from the paper;
(b) for PINNED strings (one fixed setting), the per-rectangle terms A_i of Eq. (20)
equal the coefficient of s^i in prod_j (R_j + s*C_j) (feasible DP, O(t^2)),
and sum_i A_i = 2^t [x=v][y=u]. So the true value of a level-1 subcall on a pinned
input has a poly-time closed form, although its support in the (omega,theta)
strings has C(t,i) elements per rectangle."""
import itertools, random
def bits(a): return [(a>>k)&1 for k in range(3)]
def phi(p,q):
    s=sum(p[k]*q[kk] for k in range(3) for kk in range(k,3))
    s+=q[0]*p[1]*p[2]+p[0]*q[1]*p[2]+p[0]*p[1]*q[2]
    return s
def sig(a,b): return -1 if phi(bits(a),bits(b))%2 else 1
ok=all(2*(x==v)*(y==u)==(1 if x^y^u^v==0 else 0)*(sig(x,y)*sig(v,u)+sig(x,u)*sig(v,y))
       for x,y,u,v in itertools.product(range(8),repeat=4))
print("eight-square identity (12) holds on all 4096 settings:",ok)
def Rt(x,y,u,v): return (x^y^u^v==0)*sig(x,y)*sig(v,u)
def Ct(x,y,u,v): return (x^y^u^v==0)*sig(x,u)*sig(v,y)
rng=random.Random(1); bad=0
for trial in range(2000):
    t=rng.randint(1,7)
    blk=[tuple(rng.randrange(8) for _ in range(4)) for _ in range(t)]
    if rng.random()<0.5:  # force consistency sometimes
        blk=[(x,y,y,x) for (x,y,_,_) in blk]
    R=[Rt(*b) for b in blk]; C=[Ct(*b) for b in blk]
    brute=[0]*(t+1)
    for T in itertools.product((0,1),repeat=t):
        p=1
        for j in range(t): p*=C[j] if T[j] else R[j]
        brute[sum(T)]+=p
    dp=[1]
    for j in range(t):
        new=[0]*(len(dp)+1)
        for i,c in enumerate(dp): new[i]+=c*R[j]; new[i+1]+=c*C[j]
        dp=new
    target=2**t*all(x==v and y==u for x,y,u,v in blk)
    if dp!=brute or sum(brute)!=target: bad+=1
print("A_i = DP coefficient and sum_i A_i = 2^t[x=v][y=u]; failures:",bad,"/2000")
