import random, sys
from relax import *
algs = {'BLP': blp, 'AIP_Z': aip, 'AIP_GF2': aip_gf2, 'BLP+AIP': blp_aip}
rng = random.Random(11)
for name, d in [('2-SAT', 1.8), ('1-in-3', 0.9), ('Horn', 2.0), ('XORSAT', 1.2)]:
    T = TEMPLATES[name]; st = {a: [0, 0] for a in algs}; nsat = 0; N = 200
    for _ in range(N):
        n = rng.randint(4, 9); m = max(2, int(round(d * n * rng.uniform(0.7, 1.3))))
        cons = rand_instance(T, n, m, rng); tr = brute(n, cons); nsat += tr
        for a, f in algs.items():
            r = f(n, cons); st[a][0 if r == tr else 1] += 1
    print(name, f'd={d} sat={nsat}/{N}', ' '.join(f'{a}:{c}/{N}(FA={fa})' for a, (c, fa) in st.items()), flush=True)
