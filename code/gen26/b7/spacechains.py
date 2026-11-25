import math, numpy as np, b7lib as B, chain_b7 as C
def run(label, S):
    best=None
    for name,x,e in S:
        for zname,zf in (('new',B.zeta_new),):
            cands=[C.chain_greedy(x,e,zf)]+[C.chain_param(x,e,zf,a0,b0) for a0 in np.arange(-4,16.01,1.0) for b0 in (0,.5,1,1.5,2)]
            for ch in cands:
                mu,s2=C.summarize(ch)
                if best is None or mu>best[0]: best=(mu,s2,name)
    mu,s2,name=best
    print(f"{label:46s} mu={mu:.4f} ({name}) single c1(a=2)={s2:.3f} no-saw a=2 {2*math.sqrt(2/mu):.3f} a=4 {2*math.sqrt(4/mu):.3f}",flush=True)
run("poly space: clique D&C base + R2", C.starts('clique-DC'))
run("poly space: row D&C level1 (zs=2 zs_R3 [H]) + R2", C.starts('row-DC'))
for M in (1/64, 1/16, 0.1, 0.15, 0.2, 0.25):
    S=[(f"x0={x0}",float(x0),(1-M)/3 + 1.0/2**x0) for x0 in range(4,16)]
    run(f"clique base, memory M={M:.4f}: e0=(1-M)/3+1/K", S)
    S=[(f"g=2^{lg}",2*lg,(C.RHO+2.0**-lg)*(1-M)/3+(1+ (0 if M==0.25 else 1))*B.zs_R3(lg)) for lg in np.arange(6,40.01,.5)]
    run(f"row level1, memory M={M:.4f} (rel. to U) [H]", S)
