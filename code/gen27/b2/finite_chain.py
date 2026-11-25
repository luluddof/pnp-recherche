# Tensor powers of a finite chain of realizable pairs at size N0: pieces A_s = R_s - R_{s-1}, ranks <= (p_s+p_{s-1}, q_{s-1}+q_s).
# Grouping by type costs entropy => h(mu) = log_{N0^2} sum_s r_s^mu c_s^(1-mu)  (log-sum-exp). rho = h(1/2).
import numpy as np
def analyse(N0, pairs, label, best_case=False):
    pairs=sorted(set(pairs)); # increasing p
    rs=[];cs=[]
    for (p0,q0),(p1,q1) in zip(pairs[:-1],pairs[1:]):
        if best_case: rs.append(max(p1,1)); cs.append(max(q0,1))       # optimistic: strongly nested chain
        else: rs.append(p0+p1); cs.append(q0+q1)
    rs=np.array(rs,float); cs=np.array(cs,float); L=np.log(N0**2)
    h=lambda mu: np.log(np.sum(rs**mu*cs**(1-mu)))/L
    rho=h(0.5)
    best=0
    for mu in np.linspace(0.5,0.999,500):
        hm=h(mu); lo,hi=0,0.999999
        f=lambda a:(mu**(1/(1-a))+(1-mu)**(1/(1-a)))**(1-a)
        if f(0)<hm: continue
        for _ in range(60):
            m=(lo+hi)/2
            if f(m)>=hm: lo=m
            else: hi=m
        best=max(best,lo)
    print(f"{label:38s} N0={N0} pieces={len(rs)} rho={rho:.4f} (alpha_bal={max(0,np.log2(1/rho)):.3f}) alpha*={best:.3f} h(1)={h(1.0):.3f}")
N7=[(49,0),(21,1),(20,3),(17,4),(14,5),(12,6),(11,7),(8,8)]
sym=lambda P: P+[(q,p) for p,q in P]
analyse(7, sym(N7), "N=7 exact frontier, chain (sum ranks)")
analyse(7, sym(N7), "N=7 exact frontier, strongly nested", True)
analyse(8, sym([(64,0),(28,1),(22,5),(20,6),(8,8)]), "N=8 known points (sum ranks)")
analyse(8, sym([(64,0),(28,1),(22,5),(20,6),(8,8)]), "N=8 known points, strongly nested", True)
analyse(8, sym([(64,0),(8,8)]), "N=8 octonion only (Lemma 2 check)")
