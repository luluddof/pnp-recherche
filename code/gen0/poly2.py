import itertools
B=[0,1]
def clause(signs):  # signes: 1 = positif, 0 = negatif
    return [t for t in itertools.product(B,repeat=len(signs)) if any(t[i]==s for i,s in enumerate(signs))]
L = {
 "3-SAT (8 types de clauses)":  [clause(s) for s in itertools.product(B,repeat=3)],
 "2-SAT (4 types de clauses)":  [clause(s) for s in itertools.product(B,repeat=2)],
 "XORSAT (x^y^z = 0 ou 1)":     [[t for t in itertools.product(B,repeat=3) if sum(t)%2==r] for r in B],
 "Horn-SAT":                    [clause(s) for s in [(0,0,1),(0,0,0),(1,),(0,)]],
 "1-in-3 SAT":                  [[t for t in itertools.product(B,repeat=3) if sum(t)==1]],
}
ops=list(itertools.product(B,repeat=8)); f=lambda tab,a,b,c: tab[a*4+b*2+c]
proj={(0,0,0,0,1,1,1,1),(0,0,1,1,0,0,1,1),(0,1,0,1,0,1,0,1)}
names={(0,0,0,1,0,1,1,1):"majorite",(0,1,1,0,1,0,0,1):"Maltsev x^y^z",(0,0,0,0,0,0,0,1):"ET",(0,1,1,1,1,1,1,1):"OU",(0,0,0,0,0,0,0,0):"constante 0",(1,1,1,1,1,1,1,1):"constante 1"}
def pres(tab,rel):
    rs=set(rel); k=len(rel[0])
    return all(tuple(f(tab,a[i],b[i],c[i]) for i in range(k)) in rs for a in rel for b in rel for c in rel)
for name,rels in L.items():
    P=[t for t in ops if all(pres(t,r) for r in rels) and t not in proj]
    print(f"{name:28} polymorphismes non triviaux={len(P):2}  ->  {', '.join(names.get(t,'autre') for t in P) or 'AUCUN'}")
