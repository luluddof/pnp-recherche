import numpy as np, time, sys
from starlib import *
n = tuple(int(a) for a in sys.argv[1].split(',')); rs = [int(a) for a in sys.argv[2].split(',')]; starts = int(sys.argv[3]); maxit = int(sys.argv[4])
tl = float(sys.argv[5]) if len(sys.argv) > 5 else 1e9
inst = Inst(*n); print('n', n, 'd', inst.d, 'D', inst.D, 'npar', inst.npar, 'product r', inst.prod_r, flush=True)
T0 = time.time()
for r in rs:
    fs = []; t0 = time.time(); best = (9e9, None)
    for seed in range(starts):
        x0 = inst.random_x(np.random.default_rng(1000*r + seed), [0.1, 0.03, 0.3][seed % 3])
        x, f, h = lm_solve(inst, r, x0, maxit=maxit, tlim=tl)
        fs.append(f)
        if f < best[0]: best = (f, x)
        if time.time() - T0 > 540: break
    nm = norms(inst, best[1], r); fs = np.array(fs)
    print(f"n={n} r={r}: succ(<1e-9) {int((fs<1e-9).sum())}/{len(fs)} best={fs.min():.3e} median={np.median(fs):.3e} worst={fs.max():.3e} "
          f"t={time.time()-t0:.0f}s | best: stars {np.round(nm['star_norms'],2)} |B|={nm['B_norm']:.2f} |A|={nm['A_norm']:.2f} sig_r={nm['sig_r']} sig_r+1={nm['sig_r1']:.1e}", flush=True)
    np.savez(f"sol_n{''.join(map(str,n))}_r{r}.npz", x=best[1], f=best[0], fs=fs)
