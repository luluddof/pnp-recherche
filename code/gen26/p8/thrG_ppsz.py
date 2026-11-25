import sys, math, time
sys.path.insert(0,'.')
import run_p8 as R
w=R.set_alpha(4.26); t0=time.time(); tau=0.386
lo,hi=6.06e6/1.04,6.06e6*1.04
while hi/lo>1.002:
    mid=math.sqrt(lo*hi); e,D=R.Ebest(mid,w)
    print(f"  lam={mid:.5g} E/n={e:.5f} D={D} [{time.time()-t0:.0f}s]",flush=True)
    if e<tau: hi=mid
    else: lo=mid
print(f"alpha=4.26 row G vs PPSZ: lam* in ({lo:.4g},{hi:.4g}]",flush=True)
