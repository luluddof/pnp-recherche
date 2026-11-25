import sys, math, time
sys.path.insert(0,'.')
import run_p8 as R
alpha=float(sys.argv[1]); w=R.set_alpha(alpha); t0=time.time()
for lam in (1e5,1e6,1e7,1e8,1e9,1e10):
    eg,dg=R.Ebest(lam,w,Ds=(2,3,4)); eh,dh=R.Ebest(lam,w,Ds=(2,3,4),smith=True)
    print(f"alpha={alpha} lam={lam:.0e}: G {eg*math.sqrt(lam)/w:.2f} (D={dg})  H {eh*math.sqrt(lam)/w:.2f} (D={dh})  ratio G/H {eg/eh:.4f} [{time.time()-t0:.0f}s]",flush=True)
