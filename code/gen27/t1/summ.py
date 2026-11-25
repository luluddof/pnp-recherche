import sys, pickle, glob; sys.dont_write_bytecode = True
for fn in sys.argv[1:]:
    d = pickle.load(open(fn, 'rb')); h = d['h']
    ks = [d['res'][W][1] if W in d['res'] else None for W in range(4*h)]
    wp = [d['res'][W][0] if W in d['res'] else None for W in range(4*h)]
    print(f"h={h} D={d['D']:>2} r_all={d['r_all']:>5}/{d['n']**3}: k(W)={ks}  minW'={wp}")
