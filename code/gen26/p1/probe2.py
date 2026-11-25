import sys, pickle, time; sys.dont_write_bytecode = True
from feas import feasible
h = int(sys.argv[1]); fns = sys.argv[2].split(','); Ws = [int(a) for a in sys.argv[3].split(',')]; ks = [int(a) for a in sys.argv[4].split(',')]
ds = [pickle.load(open(f, 'rb')) for f in fns]
for d in ds: d['h'] = h
print(f"h={h} sectors {fns} r_all {[d['r_all'] for d in ds]}", flush=True)
for W in Ws:
    for k in ks:
        Wp = 4*h - 1 - W - k
        if Wp < 0: continue
        t = time.time(); fs = [feasible(d, W, Wp) for d in ds]
        print(f"  W={W:>3} W'={Wp:>3} k={k}: feasible={all(f[0] for f in fs)}  deficiency per sector {[f[2]-f[1] for f in fs]}  [{time.time()-t:.0f}s]", flush=True)
