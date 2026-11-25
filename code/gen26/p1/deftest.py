import sys, pickle; sys.dont_write_bytecode = True
from feas import feasible
d = pickle.load(open(sys.argv[1], 'rb')); h = d['h']; Ws = [int(a) for a in sys.argv[2].split(',')]; kmax = int(sys.argv[3])
for W in Ws:
    out = []
    for k in range(0, kmax + 1):
        Wp = 4*h - 1 - W - k
        if Wp < 0: break
        f = feasible(d, W, Wp); out.append(f"k{k}:{f[2]-f[1]}")
    print(f"h={h} W={W}: deficiency (number of missing l's) by k: " + " ".join(out), flush=True)
