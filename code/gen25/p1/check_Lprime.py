# Brute-force check of Lemma L' (subformula line lemma) on random small instances
# Phi = phi_t (at-most-one per clique) AND M (random monotone clauses), random partial assignment rho.
# Claim: every clause (~u v ~v) with u,v free in the same clique and u live (some solution of Phi|rho has u=1)
# is NECESSARY in Phi|rho (removing it alone changes the solution set), hence lies in every equivalent subformula.
import itertools, random, sys
random.seed(int(sys.argv[1]) if len(sys.argv)>1 else 1)
def sat_set(clauses, n, rho):
    out=[]
    for x in itertools.product([0,1],repeat=n):
        if any(x[v]!=b for v,b in rho.items()): continue
        if all(any((x[abs(l)-1]==1) if l>0 else (x[abs(l)-1]==0) for l in C) for C in clauses): out.append(x)
    return set(out)
bad=0; tested=0
for trial in range(300):
    N=random.choice([2,3]); t=random.choice([3,4]); n=N*t
    cl=[]
    Q=[list(range(c*t,(c+1)*t)) for c in range(N)]
    for q in Q:
        for u,v in itertools.combinations(q,2): cl.append((-(u+1),-(v+1)))
    M=[tuple(random.sample(range(1,n+1),random.choice([1,2,3]))) for _ in range(random.randint(0,4))]
    allc=cl+M
    rho={v:random.randint(0,1) for v in random.sample(range(n),random.randint(0,3))}
    S0=sat_set(allc,n,rho)
    if not S0: continue
    for q in Q:
        free=[v for v in q if v not in rho]
        live=[u for u in free if any(s[u]==1 for s in S0)]
        for u in live:
            for v in free:
                if v==u: continue
                C=(-(min(u,v)+1),-(max(u,v)+1))
                rest=[D for D in allc if D!=C]
                tested+=1
                if sat_set(rest,n,rho)==S0: bad+=1; print("counterexample",N,t,M,rho,u,v)
print("tested",tested,"bad",bad)
