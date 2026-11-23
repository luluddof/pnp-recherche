import itertools
# Relations booleennes (ensembles de tuples autorises)
R = {
 "3-SAT  (x|y|z)":        [t for t in itertools.product([0,1],repeat=3) if any(t)],
 "XOR    (x^y^z=1)":      [t for t in itertools.product([0,1],repeat=3) if sum(t)%2==1],
 "2-SAT  (x|y)":          [t for t in itertools.product([0,1],repeat=2) if any(t)],
 "Horn   (~x|~y|z)":      [t for t in itertools.product([0,1],repeat=3) if (not t[0]) or (not t[1]) or t[2]],
 "1-in-3 (x+y+z=1)":      [t for t in itertools.product([0,1],repeat=3) if sum(t)==1],
}
# toutes les operations ternaires f:{0,1}^3->{0,1} (256), codees par leur table
ops = list(itertools.product([0,1],repeat=8))
def f(tab,a,b,c): return tab[a*4+b*2+c]
proj = {(0,0,0,0,1,1,1,1),(0,0,1,1,0,0,1,1),(0,1,0,1,0,1,0,1)}
names = {(0,0,0,1,0,1,1,1):"majorite", (0,1,1,0,1,0,0,1):"Maltsev x^y^z", (0,0,0,0,0,0,0,1):"ET", (0,1,1,1,1,1,1,1):"OU"}
def preserves(tab, rel):
    rs=set(rel)
    for t1 in rel:
        for t2 in rel:
            for t3 in rel:
                if tuple(f(tab,t1[i],t2[i],t3[i]) for i in range(len(t1))) not in rs: return False
    return True
for name,rel in R.items():
    P=[tab for tab in ops if preserves(tab,rel)]
    nontriv=[tab for tab in P if tab not in proj]
    idem=[tab for tab in nontriv if tab[0]==0 and tab[7]==1]
    known=[names[t] for t in nontriv if t in names]
    print(f"{name:22} polymorphismes ternaires={len(P):3}  non triviaux idempotents={len(idem):3}  connus: {', '.join(known) or 'aucun'}")
