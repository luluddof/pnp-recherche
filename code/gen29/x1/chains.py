# Finite-chain support functions (B2 formula): tensor powers of a chain of realizable pairs at base dimension D0
# (D0 = full rank of the operator in each direction): pieces r_s = p_{s-1}+p_s, c_s = q_{s-1}+q_s (sum) or
# optimistic r_s = p_s, c_s = q_{s-1}; h(mu) = log_{D0} sum_s r_s^mu c_s^(1-mu). alpha* over mu in (0,1) (both corners).
import numpy as np
def h_fun(D0,pairs,opt):
    pairs=sorted(set(pairs)); rs=[];cs=[]
    for (p0,q0),(p1,q1) in zip(pairs[:-1],pairs[1:]):
        if opt: rs.append(max(p1,1)); cs.append(max(q0,1))
        else: rs.append(p0+p1); cs.append(q0+q1)
    rs=np.array(rs,float); cs=np.array(cs,float); L=np.log(D0)
    return lambda mu: np.log(np.sum(rs**mu*cs**(1-mu)))/L
def alpha_star(h):
    best=0
    for mu in np.concatenate([np.linspace(0.001,0.999,999)]):
        hm=h(mu); lo,hi=0,0.999999
        f=lambda a:(mu**(1/(1-a))+(1-mu)**(1/(1-a)))**(1-a)
        if f(0)<hm: continue
        for _ in range(60):
            m=(lo+hi)/2
            if f(m)>=hm: lo=m
            else: hi=m
        best=max(best,lo)
    return best
data={
 'SWAP N=2':(4,[(4,0),(1,1),(0,4)]),
 'SWAP N=3':(9,[(9,0),(3,1),(1,3),(0,9)]),
 'SWAP N=4':(16,[(16,0),(6,1),(5,3),(4,4),(3,5),(1,6),(0,16)]),
 'SWAP N=7 (exact frontier)':(49,[(49,0),(21,1),(20,3),(17,4),(14,5),(12,6),(11,7),(8,8),(7,11),(6,12),(5,14),(4,17),(3,20),(1,21),(0,49)]),
 'X_{2,3}':(6,[(6,0),(3,1),(2,2),(1,3),(0,6)]),
 'X_{2,4} (=3-cycle 2-vs-1, N=2)':(8,[(8,0),(5,1),(2,2),(1,5),(0,8)]),
 'X_{3,4}':(12,[(12,0),(6,1),(5,2),(4,3),(3,4),(2,5),(1,6),(0,12)]),
 'X_{2,5}':(10,[(10,0),(7,1),(4,2),(3,3),(2,4),(1,7),(0,10)]),
 'P_3 ring N=2':(8,[(8,0),(4,1),(2,4),(0,8)]),
}
import sys
extra=eval(sys.argv[1]) if len(sys.argv)>1 else {}
data.update(extra)
for k,(D0,pr) in data.items():
    for opt in (False,True):
        h=h_fun(D0,pr,opt)
        print(f"{k:34s} {'optimistic' if opt else 'sum-ranks '} rho=h(1/2)={h(0.5):.4f} alpha_bal={max(0,np.log2(1/h(0.5))):.3f} alpha*={alpha_star(h):.3f}")
