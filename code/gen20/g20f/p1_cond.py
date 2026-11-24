# (3.2): every partition has M_i >= phi * M * min(1, 2(1-beta)), phi = min_{0<vol S<=M} cut(S)/vol(S), beta = max vol(g)/2M
import numpy as np, random, itertools
random.seed(3); rng=np.random.default_rng(3)
viol=0; tests=0; tight=1e9; tests_sep=0; viol_sep=0
for trial in range(300):
    n=random.randint(3,13); M=random.randint(n,3*n)
    wires=[(random.randrange(n),random.randrange(n)) for _ in range(M)]
    deg=np.zeros(n,int)
    for h,t in wires: deg[h]+=1; deg[t]+=1
    masks=np.arange(1,2**n-1); bits=((masks[:,None]>>np.arange(n))&1).astype(bool)
    vol=bits@deg
    cut=np.zeros(len(masks),int)
    for h,t in wires: cut+=bits[:,h]!=bits[:,t]
    ok=(vol>0)&(vol<=M)
    if not ok.any(): continue
    phi=(cut[ok]/vol[ok]).min()
    for _ in range(40):
        r=random.randint(1,6); lab=np.array([min(r-1,int(rng.exponential(random.choice([0.3,1,3])))) for _ in range(n)])
        Mi=sum(lab[h]!=lab[t] for h,t in wires)
        beta=max(deg[lab==g].sum() for g in range(r))/(2*M)
        bound=phi*M*min(1,2*(1-beta)); tests+=1
        if Mi < bound-1e-9: viol+=1
        if bound>0: tight=min(tight,Mi/bound)
        # separator claim: merging groups -> bipartition, both sides vol >= min(1-beta,1/3)*2M (largest-first greedy)
        vg=sorted([deg[lab==g].sum() for g in range(r)],reverse=True)
        if vg[0]>=2*M/3: side=vg[0]
        else:
            side=0
            for v in vg:
                if side<2*M/3 and side+v<=2*M*2/3+1e-9 or side< 2*M/3: side+=v
                if side>=2*M/3-1e-9 or side>=2*M/3: break
        # simple existence check by brute force over group subsets
        best=max(min(sum(vg[j] for j in range(r) if (m>>j)&1), 2*M-sum(vg[j] for j in range(r) if (m>>j)&1)) for m in range(2**r))
        tests_sep+=1
        if best < min(1-beta,1/3)*2*M-1e-9: viol_sep+=1
print('(3.2) tests',tests,'violations',viol,' min M_i/bound =',round(tight,4))
print('(3.3) balanced-bipartition-from-groups tests',tests_sep,'violations',viol_sep)
