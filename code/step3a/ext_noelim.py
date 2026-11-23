"""Same as ext.py but CaDiCaL with elim=0 (no bounded variable elimination),
so extension variables cannot be eliminated by the solver's preprocessing."""
import sys, json, time
from multiprocessing import Pool
from pysat.solvers import Solver
import ext

def solve(n, cls, proof=False):
    with Solver(name='cadical195', with_proof=proof) as s:
        s.solver.configure({'elim': 0})
        s.append_formula(cls)
        r = s.solve(); st = s.accum_stats()
        pl = len(s.get_proof()) if (proof and not r) else None
    return r, st['conflicts'], pl

ext.solve = solve

if __name__ == '__main__':
    ns = [int(x) for x in sys.argv[1].split(',')]
    per = int(sys.argv[2]); out = sys.argv[3]
    with Pool(4) as pool, open(out, 'a') as f:
        for n in ns:
            got, seed = 0, 0
            while got < per:
                batch = [(n, s) for s in range(seed, seed + 2 * per)]
                seed += 2 * per
                for res in pool.imap(ext.job, batch):
                    if res and got < per:
                        got += 1; f.write(json.dumps(res) + '\n'); f.flush()
            print('done n=', n, flush=True)
