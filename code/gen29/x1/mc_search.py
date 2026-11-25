# Exact search (integer entries in {-1,0,1}) for decompositions T = R + C with C in a column format,
# minimizing rank_rows(T - C). Annealing on single entries; final check with flint (rank over Z, mod 2 primes).
# usage: python3 mc_search.py rect N1 N2 q target_p [secs]   |  cyc3 N b12,b23,b31 target_p [secs]  | copy N q target_p
import numpy as np, sys, time, json, flint
from vp import Rect, Cyc3, Copy

def rank_np(M): return np.linalg.matrix_rank(M.astype(float), tol=1e-8)

def exact_ranks(M):
    Z = flint.fmpz_mat([[int(v) for v in r] for r in M]); out = {'Q': Z.rank()}
    for pr in (1000003, 998244353): out[pr] = flint.nmod_mat([[int(v) % pr for v in r] for r in M], pr).rank()
    return out

def run(prob, target, secs, seed, density=0.3):
    rng = np.random.default_rng(seed)
    z = rng.choice([-1, 0, 1], size=prob.n, p=[density / 2, 1 - density, density / 2]).astype(float)
    cur = rank_np(prob.rows(prob.T - prob.C(z)[0])); best = (cur, z.copy()); T0 = 1.0; t0 = time.time(); it = 0
    while time.time() - t0 < secs and best[0] > target:
        it += 1; temp = T0 * max(0.02, 1 - (time.time() - t0) / secs)
        z2 = z.copy(); k = rng.integers(prob.n, size=rng.integers(1, 3))
        z2[k] = rng.choice([-1, 0, 1], size=len(k))
        r2 = rank_np(prob.rows(prob.T - prob.C(z2)[0]))
        if r2 <= cur or rng.random() < np.exp(-(r2 - cur) / temp):
            z, cur = z2, r2
            if cur < best[0]: best = (cur, z.copy())
    return best

if __name__ == '__main__':
    kind = sys.argv[1]
    if kind == 'rect':
        N1, N2, q, tp = map(int, sys.argv[2:6]); prob = Rect(N1, N2, q); rest = sys.argv[6:]; tag = f"rect_{N1}x{N2}_q{q}"
    elif kind == 'cyc3':
        N = int(sys.argv[2]); b = tuple(map(int, sys.argv[3].split(','))); tp = int(sys.argv[4]); prob = Cyc3(N, b)
        rest = sys.argv[5:]; tag = f"cyc3_N{N}_b{b[0]}{b[1]}{b[2]}"
    elif kind == 'copy':
        N, q, tp = map(int, sys.argv[2:5]); prob = Copy(N, q); rest = sys.argv[5:]; tag = f"copy_N{N}_q{q}"
    secs = float(rest[0]) if rest else 60
    best = (10 ** 9, None); seed = 0; t0 = time.time()
    while time.time() - t0 < secs and best[0] > tp:
        b_ = run(prob, tp, min(20, secs), seed); seed += 1
        if b_[0] < best[0]: best = b_
    z = best[1]; M = prob.rows(prob.T - prob.C(z)[0]).astype(int)
    er = exact_ranks(M)
    # column-format check: for rect/copy the OSR of C exactly; for cyc3 the bonds are structural
    extra = ''
    if kind in ('rect', 'copy'):
        Cz = prob.C(z)[0]
        if kind == 'rect': Ccol = Cz.transpose(0, 2, 1, 3).reshape(prob.N1 * prob.N2, prob.N2 * prob.N1)
        else: Ccol = Cz.transpose(0, 2, 4, 1, 3).reshape(prob.N ** 3, prob.N ** 2)
        extra = f" exact OSR(C) = {exact_ranks(Ccol.astype(int))}"
    print(f"{tag}: best rank found {best[0]} (target {tp}); exact rank_rows(T-C) = {er}{extra}", flush=True)
    if best[0] <= tp:
        json.dump({'kind': kind, 'args': sys.argv[2:], 'z': [int(v) for v in z], 'ranks': {str(k): v for k, v in er.items()}},
                  open(f"exact_{tag}_p{er['Q']}.json", 'w'))
