"""one start: run1.py n01,n02,n12 r seed core_scale maxit tlim [init.npz]   -> appends a line to runs.log and saves ck_*.npz"""
import numpy as np, time, sys, json
from starlib import *
n = tuple(int(a) for a in sys.argv[1].split(',')); r = int(sys.argv[2]); seed = int(sys.argv[3]); cs = float(sys.argv[4])
maxit = int(sys.argv[5]); tl = float(sys.argv[6]); init = sys.argv[7] if len(sys.argv) > 7 else None
cgm = int(sys.argv[8]) if len(sys.argv) > 8 else 40
inst = Inst(*n); tag = f"n{''.join(map(str,n))}_r{r}_s{seed}" + ("_cont" if init else "")
if init:
    x0 = np.load(init)['x']
    if cs > 0: x0 = x0 + cs * np.random.default_rng(seed).normal(size=x0.size) * np.linalg.norm(x0) / np.sqrt(x0.size)
else: x0 = inst.random_x(np.random.default_rng(seed), cs)
t0 = time.time()
x, f, h = lm_solve(inst, r, x0, maxit=maxit, tlim=tl, verbose=True, cg_max=cgm)
nm = norms(inst, x, r)
np.savez(f"ck_{tag}.npz", x=x, f=f, hist=np.array(h))
line = dict(n=n, r=r, seed=seed, cs=cs, init=init, f=f, its=len(h)-1, time=round(time.time()-t0,1), **nm)
open('runs.log', 'a').write(json.dumps(line) + "\n"); print(json.dumps(line))
