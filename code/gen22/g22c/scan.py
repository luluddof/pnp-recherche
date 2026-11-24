# generic scan: for each (family, N) compute p at the best scalar (best_scale; fallback scales 1, 2, 1/2) and record q.
import sys, time, csv; sys.path.insert(0, '.')
import numpy as np, scipy.sparse as sp
from ev import *
from fams import *
from cd import table
def evaluate(N, terms, try_best=True):
    res = []
    for num, den in [(1, 1), (2, 1), (1, 2), (-1, 1)]:
        t2 = [(G, Gp * den) for G, Gp in terms]
        res.append((rank(build(N, t2, d=num))[0], f'{num}/{den}'))
    if try_best and N * N <= 6000:
        try:
            bs = best_scale(N, terms)
            if bs:
                for r, lam, _ in bs:
                    if r is not None: res.append((r, f'{lam[0]}/{lam[1]}'))
        except Exception as e: print('best_scale fail', e)
    return min(res)
def instances(name, N):
    """yield (label, terms) for admissible N"""
    if name == 'cyc_struct':
        T, inv = grp_cyclic(N); yield 'f=1', struct_slices(T, np.ones((N, N), int))
        rng = np.random.default_rng(N); f = rng.choice([-1, 1], (N, N)); yield 'f=rand', struct_slices(T, f)
    if name == 'elab2_struct' and N & (N - 1) == 0 and N > 1:
        n = N.bit_length() - 1; T, inv = grp_elab(2, n)
        S, K = table(n); yield 'CD', struct_slices(T, S)
        # Clifford cocycle (-1)^{sum_{i>j} a_i c_j}
        el = np.array([[(v >> k) & 1 for k in range(n)] for v in range(N)])
        B = np.tril(np.ones((n, n), int), -1); f = (-1) ** ((el @ B @ el.T) % 2); yield 'clifford', struct_slices(T, f)
    if name == 'elab3_struct':
        n = round(np.log(N) / np.log(3))
        if 3 ** n == N and n >= 1:
            T, inv = grp_elab(3, n); yield 'f=1', struct_slices(T, np.ones((N, N), int))
    if name == 'dihedral_struct' and N % 2 == 0 and N >= 4:
        T, inv = grp_dihedral(N // 2); yield 'f=1', struct_slices(T, np.ones((N, N), int))
    if name == 'heis_struct':
        for pr in (2, 3, 5):
            if pr ** 3 == N:
                T, inv = grp_heis(pr); yield 'f=1', struct_slices(T, np.ones((N, N), int))
    if name == 'translations':
        T, inv = grp_cyclic(N)
        for mode in ('same', 'T', 'R'):
            yield 'cyc-' + mode, translations(T, inv, range(N), mode)
        if N & (N - 1) == 0 and N > 1:
            T2, i2 = grp_elab(2, N.bit_length() - 1)
            for mode in ('same',): yield 'elab2-' + mode, translations(T2, i2, range(N), mode)
        if N % 2 == 0 and N >= 4:
            T3, i3 = grp_dihedral(N // 2)
            for mode in ('same', 'T', 'R'): yield 'dih-' + mode, translations(T3, i3, range(N), mode)
    if name == 'field' and isprime(N) and N > 2:
        r = N; c = chi(r); sq = [a for a in range(1, r) if c[a] == 1]
        for A, lab in ((list(range(1, r)), 'all'), (sq, 'sq')):
            for mode in ('T', 'same', 'inv'):
                for sm in ('none', 'chi'):
                    yield f'mult-{lab}-{mode}-{sm}', mult_family(r, A, mode, sm)
    if name == 'affine' and isprime(N) and N > 2 and N <= 31:
        r = N; c = chi(r); sq = [a for a in range(1, r) if c[a] == 1]
        for mode in ('T', 'same'):
            for sm in ('none', 'chi'):
                yield f'aff-sq-{mode}-{sm}', affine_family(r, [(a, b) for a in sq for b in range(r)], mode, sm)
if __name__ == '__main__':
    name = sys.argv[1]; Ns = list(map(int, sys.argv[2].split(','))); tb = len(sys.argv) < 4
    w = csv.writer(open(f'scan_{name}.csv', 'a'))
    for N in Ns:
        for lab, terms in instances(name, N):
            t = time.time(); p, sc = evaluate(N, terms, tb)
            q = len(terms); print(name, lab, N, 'q', q, 'p', p, 'scale', sc, 'N^2', N * N, '%.1fs' % (time.time() - t), flush=True)
            w.writerow([name, lab, N, q, p, sc])
