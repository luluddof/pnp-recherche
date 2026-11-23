"""Experiment runner. usage: python3 run.py <family>  -> prints table rows."""
import sys, random, time
from common import *
import candA_pc as A
from candBCD import *

fam = sys.argv[1]
rng = random.Random(12345)


def instances():
    if fam == 'rand':
        for n in (30, 50, 75, 100, 150, 200):
            for _ in range(20 if n <= 100 else 10):
                yield n, rand3sat(n, round(4.26 * n), rng)
    elif fam == 'tseitin':
        for nv in (10, 20, 30, 40, 60):
            for odd in (True, False):
                for _ in range(4):
                    yield nv, tseitin(nv, rng, odd)
    elif fam == 'php':
        for h in (3, 4, 5, 6, 7):
            yield h, php(h)
    elif fam == 'xor':
        for n in (30, 60, 90):
            for ratio in (0.6, 1.0):
                for _ in range(5):
                    yield n, rand3xor(n, round(ratio * n), rng)


from collections import defaultdict
stats = defaultdict(lambda: [0, 0, 0, 0])  # (size,cand) -> [sat_ok, sat_tot, unsat_ok, unsat_tot]
for size, (n, cls) in instances():
    t = truth(n, cls)
    res = {}
    # A: degree-3 PC decimation (full decision procedure); too slow beyond ~60 vars in Python
    if n <= 60:
        r, a = A.decide(n, cls, 3)
        res['A:PC3'] = (r == 'SAT' and check(cls, a)) if t else (r == 'UNSAT')
    # A': degree-4 PC refutation only (small n)
    if not t and n <= 30:
        res['A:PC4-ref'] = A.refute(n, cls, 4)
    # B: LP rounding + WalkSAT(20 n^2 flips) for SAT; cluster-LP infeasibility for UNSAT
    if t:
        a = candB(n, cls, rng); res['B:LP+WS'] = a is not None and check(cls, a)
    elif n <= 30:
        res['B:clusterLP'] = not cluster_lp_feasible(n, cls)
    else:
        res['B:basicLP'] = basic_lp(n, cls, rng) is None
    # C: BP-guided decimation (cannot certify UNSAT)
    if t:
        a = bp_decimate(n, cls, rng); res['C:BPdec'] = a is not None and check(cls, a)
    else:
        res['C:BPdec'] = False
    # D: spectral refutation (cannot find solutions)
    if not t:
        if all(len(c) == 3 for c in cls):
            res['D:spectral'] = bool(spectral_refute(n, cls)[0])
    for k, ok in res.items():
        s = stats[(size, k)]
        if t: s[0] += ok; s[1] += 1
        else: s[2] += ok; s[3] += 1
    print(fam, size, n, 'SAT' if t else 'UNSAT', res, flush=True)

print('=== SUMMARY', fam)
for (size, k), s in sorted(stats.items()):
    f = lambda a, b: f'{a}/{b}' if b else '-'
    print(f'{fam:8s} size={size:4d} {k:12s} SAT {f(s[0], s[1]):7s} UNSAT {f(s[2], s[3])}')
