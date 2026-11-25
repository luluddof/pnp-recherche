import sys, pickle; sys.dont_write_bytecode = True
from feas import rank_mod, allowed
d = pickle.load(open(sys.argv[1], 'rb')); h = d['h']; P = d['P']
wset = set(int(a) for a in sys.argv[2].split(','))
def feas(W, Wp):
    cols = [v for (mu, e), v in d['cols'].items() if sum(mu) in wset and allowed(h, mu, e, W, Wp)]
    allc = [v for (mu, e), v in d['cols'].items() if sum(mu) in wset]
    return rank_mod(cols, P) == rank_mod(cols + d['rhs'], P), rank_mod(allc, P), rank_mod(allc + d['rhs'], P)
print("weights", sorted(wset), "span of all cols contains all l(x)Id:", feas(4*h, 4*h)[1:] )
for W in [h, 2*h - 1, 2*h, 3*h - 1]:
    res = []
    for k in range(0, 6):
        Wp = 4*h - 1 - W - k
        res.append(feas(W, Wp)[0])
    print(f"W={W}: feasible for k=0..5:", res, flush=True)
