import sys, math, time
sys.path.insert(0,'.')
import run_p8 as R
lam=float(sys.argv[1]); w=R.set_alpha(4.26); t0=time.time()
for D in (3,4,5):
    eh=R.optG(lam,D,w,smith=True)
    print(f"H lam={lam:.0e} D={D}: E sqrt(lam)/(km) = {eh*math.sqrt(lam)/w:.3f}  [{time.time()-t0:.0f}s]",flush=True)
