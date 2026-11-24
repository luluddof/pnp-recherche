#!/usr/bin/env python3
"""G18S r3: exact GF(2) search (p=1, kappa=1, 0/1 tables) on the small cases: try to beat the A'/A'' bound (refutation
attempt) and find the least N with a mod-2 DF identity.  Each (N, direction profile) is a separate SAT call with a wall-clock cap."""
import itertools, sys, time, multiprocessing as mp
from dflib import *
from r2_small import INST
CAP = int(sys.argv[1]) if len(sys.argv) > 1 else 120
def job(args):
    u, wires, dirs = args
    t = time.time(); r, tabs = sat_search(u, wires, dirs)
    ok = verify_gf2(u, wires, dirs, tabs) if r == 'SAT' else None
    return dirs, r, ok, round(time.time() - t, 1), tabs
def run(u, wires, dirs, cap):
    with mp.Pool(1) as pool:
        a = pool.apply_async(job, ((u, wires, dirs),))
        try: return a.get(cap)
        except mp.TimeoutError:
            pool.terminate(); return dirs, 'TIMEOUT', None, cap, None
if __name__ == "__main__":
    only = sys.argv[2:] 
    for name, u, wires in INST:
        if only and not any(o in name for o in only): continue
        b = bounds(u, wires); lb = 2 ** max(b['Ap'], b['App'])
        print(f"{name}: M={len(wires)} proved N >= {lb}", flush=True)
        found = False
        for N in range(max(1, lb - 1), 2 ** len(wires) + 1):
            res = []
            for dirs in itertools.combinations_with_replacement(range(u), N):
                d, r, ok, tt, tabs = run(u, wires, dirs, CAP); res.append((d, r, ok, tt))
                if r == 'SAT': found = True; break
            print(f"   N={N}: " + "; ".join(f"{d}:{r}{'' if ok is None else '(verified)' if ok else '(BAD)'} {tt}s" for d, r, ok, tt in res), flush=True)
            if found or any(r == 'TIMEOUT' for _, r, _, _ in res) and N >= lb + 3: break
