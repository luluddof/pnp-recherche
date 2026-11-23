"""Exact (or bounded) B2 circuit size of named 5/6-input functions + random 5-input functions."""
import sys, time, random, itertools
from synth import circuit_size
TO = float(sys.argv[1]) if len(sys.argv) > 1 else 120

def fn(n, pred):
    return sum(1 << t for t in range(1 << n) if pred([(t >> j) & 1 for j in range(n)]))

# triangle in K4: edges 01,02,03,12,13,23 -> inputs 0..5
E = {(0, 1): 0, (0, 2): 1, (0, 3): 2, (1, 2): 3, (1, 3): 4, (2, 3): 5}
tri = lambda x: any(x[E[(a, b)]] and x[E[(a, c)]] and x[E[(b, c)]] for a, b, c in itertools.combinations(range(4), 3))
# 2-SAT-like / NP-flavoured: does the 5-bit string satisfy a fixed 3-CNF?  Instead: "exists y in {0,1}^2 with
# CNF(x,y)" is too ad hoc; use clique (triangle), and a 5-input 'subset-sum' style: x0+2x1+x2+2x3+x4 == 4
named = [
    ('PAR5', 5, lambda x: sum(x) % 2),
    ('MAJ5', 5, lambda x: sum(x) >= 3),
    ('TH2_5 (>=2)', 5, lambda x: sum(x) >= 2),
    ('EXACT2_5', 5, lambda x: sum(x) == 2),
    ('MOD3_5', 5, lambda x: sum(x) % 3 == 0),
    ('SUBSETSUM5 (1,2,1,2,1 ->4)', 5, lambda x: x[0] + 2 * x[1] + x[2] + 2 * x[3] + x[4] == 4),
    ('TRIANGLE_K4 (6 edges)', 6, tri),
]
random.seed(1)
for r in range(4):
    named.append((f'RAND5_{r}', 5, (lambda f: (lambda x: (f >> sum(b << j for j, b in enumerate(x))) & 1))(random.getrandbits(32))))
which = sys.argv[2:] if len(sys.argv) > 2 else None
for name, n, p in named:
    if which and not any(w in name for w in which): continue
    t0 = time.time()
    s, exact = circuit_size(fn(n, p), n, timeout=TO, verbose=True)
    print(f'{name:30s} n={n}  size {"=" if exact else ">="} {s}   ({time.time()-t0:.0f}s)', flush=True)
