import sys, pickle, json; sys.dont_write_bytecode = True
from feas import probe
fn = sys.argv[1]; Ws = [int(a) for a in sys.argv[2].split(',')]; ks = [int(a) for a in sys.argv[3].split(',')]
Dmax = int(sys.argv[4]) if len(sys.argv) > 4 else None
d = pickle.load(open(fn, 'rb')); print(f"{fn}: h={d['h']} D={d['D']} c={d['c']} rows={d['n']**2*d['c']} Dmax={Dmax}", flush=True)
probe(d, Ws, ks, Dmax)
