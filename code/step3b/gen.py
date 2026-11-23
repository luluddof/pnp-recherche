import sys, random, os
sys.path.insert(0, '/tmp/claude-0/-home-user/2d12e629-e6d2-542b-bfc0-8b766f57383f/scratchpad/agent1')
from common import rand3sat, tseitin
from pysat.solvers import Cadical153
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'inst'); os.makedirs(D, exist_ok=True)
def issat(cls):
    s = Cadical153(bootstrap_with=cls); r = s.solve(); s.delete(); return r
def write(path, n, cls):
    with open(path, 'w') as f:
        f.write(f"p cnf {n} {len(cls)}\n")
        for c in cls: f.write(" ".join(map(str, c)) + " 0\n")
for dens in (4.26, 5.0):
    for n in (20, 30, 40, 50, 60):
        got = 0; seed = 0
        while got < 10:
            rng = random.Random(1000003 * n + int(dens * 100) * 7 + seed); seed += 1
            _, cls = rand3sat(n, round(dens * n), rng)
            if issat(cls): continue
            write(f"{D}/r{dens}_n{n}_{got}.cnf", n, cls); got += 1
# a few SAT ones for model checking
for i in range(5):
    rng = random.Random(77 + i); _, cls = rand3sat(40, 120, rng)
    assert issat(cls); write(f"{D}/sat_n40_{i}.cnf", 40, cls)
for nv in (10, 16, 20, 26, 30, 36, 42):
    for i in range(5):
        rng = random.Random(500 * nv + i); ne, cls = tseitin(nv, rng)
        assert not issat(cls); write(f"{D}/ts_v{nv}_{i}.cnf", ne, cls)
print("ok")
