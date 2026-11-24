"""G10d test 4: Thm B pipeline with the theorem's own parameters (h, l, nbar, beta), class C_t, random checks
chosen greedily until X=0, Lemma 12 truncations (M' = nu-1+r) for both X (mod 2^beta) and F (mod 2^(h+1))."""
import math, random, sys
import networkx as nx
def sigma(i): return 1 if i==0 else (-1)**i*2**(i-1)
def trunc_sum(qs, Mp):
    # sum over (Pi_1..Pi_r), Pi_j subset of a q_j-set, total size <= Mp, of prod sigma_{|Pi_j|}
    poly=[1]
    for q in qs:
        pq=[sigma(i)*math.comb(q,i) for i in range(q+1)]
        new=[0]*min(len(poly)+q, Mp+1)
        for i,a in enumerate(poly):
            if a==0: continue
            for j,b in enumerate(pq):
                if i+j<=Mp: new[i+j]+=a*b
        poly=new
    return sum(poly)
rng=random.Random(11); allok=True; runs=0; nontrivial=0
for (n,d) in [(10,3),(12,3),(12,4),(13,4),(14,3),(12,7)]:
    G=nx.random_regular_graph(d,n,seed=rng.randrange(10**6)) if (n*d)%2==0 else nx.gnm_random_graph(n,n*d//2,seed=rng.randrange(10**6))
    edges=list(G.edges()); m=len(edges)
    sets={}
    for mask in range(1<<n):
        sets.setdefault(bin(mask).count('1'),[]).append(sum(1<<i for i,(u,v) in enumerate(edges) if mask>>u&1 and mask>>v&1))
    for t in range(1,n//2+2):
        C=math.comb(n,t); h=max(1,math.ceil(math.log2(C)))
        if h>m: continue
        l=math.log2(4*m/h); nb=max(2,math.floor(h/l))
        it=sum(1 for V in sets[t] if V==0)
        beta=h+math.ceil(math.log2(math.comb(m,nb+1)))+1 if nb+1<=m else h+1
        rbar=0
        while not 3**rbar*2**beta<4**rbar: rbar+=1
        M=rbar+beta-1
        checks=[]
        def X(chk): return sum(math.comb(bin(V).count('1'),nb+1) for V in sets[t] if all(bin(V&c).count('1')%2==0 for c in chk))
        x=X(checks); x0=x
        assert x<=2**(beta-1)
        v=1
        while 2*v<=nb: v*=2
        while x>0:
            best=None
            for _ in range(12):
                c=sum(1<<e for e in range(m) if rng.random()<1.0/v)  # rate 2/v subsample, then a uniformly random parity (Lemma 13)
                xc=X(checks+[c])
                if best is None or xc<best[0]: best=(xc,c)
            checks.append(best[1]); x=best[0]
            if len(checks)>4*rbar+5: break
        r=len(checks)
        ok_r = r<=rbar
        # F truncated
        Ftr=0; Xtr=0
        for V in sets[t]:
            qs=[bin(V&c).count('1') for c in checks]
            eta=sum((-1)**j*math.comb(bin(V).count('1'),j) for j in range(nb+1))
            Ftr+=eta*trunc_sum(qs,h+r)            # nu = h+1, M' = nu-1+r
            Xtr+=math.comb(bin(V).count('1'),nb+1)*trunc_sum(qs,beta-1+r)
        okF=(Ftr%2**(h+1)==it); okX=(Xtr%2**beta==0)
        # the theorem uses M >= both
        FtrM=sum(sum((-1)**j*math.comb(bin(V).count('1'),j) for j in range(nb+1))*trunc_sum([bin(V&c).count('1') for c in checks],M) for V in sets[t])
        okM=(FtrM%2**(h+1)==it)
        runs+=1; nontrivial+= (x0>0)
        allok&=okF and okX and okM and ok_r and x==0
        print(f"n={n} d={d} m={m} t={t} h={h} nbar={nb} beta={beta} rbar={rbar} M={M} X0={x0} r={r} i_t={it} F_trunc mod 2^(h+1)={Ftr%2**(h+1)} okF={okF} okX={okX} okM={okM} r<=rbar={ok_r}")
        sys.stdout.flush()
print("runs",runs,"with X0>0:",nontrivial,"ALL OK:",allok)
