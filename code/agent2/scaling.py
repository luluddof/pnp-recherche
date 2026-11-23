"""Barrier tests: (1) largeness on random functions, (2) soundness tests against functions with KNOWN small circuits.
Known B2 upper bounds: PAR_n <= n-1, AND_n <= n-1, IP_n (n=2m) <= 2m-1, MAJ_n <= ~4.5n (O(n) via adder), so a sound
lower-bound measure mu must satisfy mu(f) <= c * (upper bound) for these."""
import numpy as np, random, math, itertools, time
from measures import anf_degree, total_influence, max_partition_rank, fourier_entropy, spectral_L1

def d2_count(tt):
    """NEW candidate 'mixed-derivative count': D2(f) = #{pairs i<j : d_i d_j f is non-constant}
    (= sum over pairs S of log2 #(subfunctions on S modulo affine functions on S)).  Range 0..C(n,2): super-linear scale."""
    t = np.asarray(tt, dtype=np.int8); n = int(math.log2(len(t))); idx = np.arange(len(t)); c = 0
    for i in range(n):
        for j in range(i + 1, n):
            d = t ^ t[idx ^ (1 << i)] ^ t[idx ^ (1 << j)] ^ t[idx ^ (1 << i) ^ (1 << j)]
            c += int(d.min() != d.max())
    return c

def d3_count(tt):
    t = np.asarray(tt, dtype=np.int8); n = int(math.log2(len(t))); idx = np.arange(len(t)); c = 0
    for i, j, k in itertools.combinations(range(n), 3):
        d = np.zeros_like(t)
        for m in range(8):
            d ^= t[idx ^ ((m & 1) << i) ^ (((m >> 1) & 1) << j) ^ (((m >> 2) & 1) << k)]
        c += int(d.min() != d.max())
    return c

def table(n, pred):
    return np.array([pred([(x >> j) & 1 for j in range(n)]) for x in range(1 << n)], dtype=np.int8)

M = dict(anf_deg=anf_degree, tot_inf=total_influence, fourier_H=fourier_entropy, log_L1=spectral_L1,
         part_logrank=max_partition_rank, D2=d2_count, D3=d3_count)
if __name__ == '__main__':
    rng = np.random.default_rng(0)
    for n in (4, 6, 8, 10):
        fams = {'PAR': (table(n, lambda x: sum(x) % 2), n - 1),
                'AND': (table(n, lambda x: int(all(x))), n - 1),
                'MAJ': (table(n, lambda x: int(sum(x) > n // 2)), None),
                'IP': (table(n, lambda x: sum(x[2 * i] & x[2 * i + 1] for i in range(n // 2)) % 2), n - 1)}
        for r in range(3):
            fams[f'RAND{r}'] = (rng.integers(0, 2, 1 << n).astype(np.int8), None)
        print(f'--- n={n}   (columns: ' + ', '.join(M) + ')')
        for name, (t, ub) in fams.items():
            vals = [M[k](t) for k in M]
            print(f'{name:6s} ub={str(ub):4s} ' + ' '.join(f'{v:7.2f}' for v in vals), flush=True)
