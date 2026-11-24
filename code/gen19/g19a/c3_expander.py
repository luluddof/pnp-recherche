#!/usr/bin/env python3
"""G19a c3: A'' (flattening bound, LP form) on an expander for (a) product labelling {0,1}^u with random balanced classes,
(b) circular labelling: 2u arcs, directions = the u diameter cuts.  r = 2 groups per direction.
A''_LP = max_{s in +-1^u} max_{x in [-1,1]^E} min_i s_i sum_e x_e delta_i(e),  delta_i(e) = (l_i(head) - l_i(tail))/2.
(the integer A'' is within u of A''_LP).  Also checks: sum_i M_i = sum_e d_e, M_i >= bisection-type bound phi*M."""
import itertools, sys, time
import numpy as np, networkx as nx
from scipy.optimize import linprog

def app_lp(delta):
    u, M = delta.shape; best = -1
    for s in itertools.product((1, -1), repeat=u - 1):
        s = np.array((1,) + s, dtype=float)
        A = np.hstack([-(delta * s[:, None]), np.ones((u, 1))])
        res = linprog(np.r_[np.zeros(M), -1.0], A_ub=A, b_ub=np.zeros(u), bounds=[(-1, 1)] * M + [(None, None)], method='highs')
        best = max(best, -res.fun)
    return best

def run(n, d, seed, us_prod, us_circ):
    rng = np.random.default_rng(seed)
    G = nx.random_regular_graph(d, n, seed=seed); E = list(G.edges()); M = len(E)
    lam = np.sort(np.linalg.eigvalsh(nx.to_numpy_array(G)))[::-1]
    phi_lb = (1 - lam[1] / d) / 2    # Cheeger: cut(S) >= phi_lb * vol(S) for vol(S) <= vol/2
    print(f"random {d}-regular graph n={n} M={M} seed={seed}: lambda_2={lam[1]:.3f}, conductance >= (1-lambda_2/d)/2 = {phi_lb:.3f}")
    def report(name, L, u):
        # L: u x n array of +-1 (balanced)
        delta = np.array([[(L[i, b] - L[i, a]) / 2 for (a, b) in E] for i in range(u)])
        Mi = np.abs(delta).sum(1); de = np.abs(delta).sum(0)
        assert abs(Mi.sum() - de.sum()) < 1e-9 and Mi.min() >= phi_lb * M - 1e-9
        t = time.time(); a = app_lp(delta)
        frac1 = (de <= 1).mean(); symf = sum(np.mean(de == k) * k * __import__('math').comb(k - 1, (k - 1) // 2) / 2 ** (k - 1) for k in range(1, u + 1)) / u
        print(f"  {name:9s} u={u:2d}: min M_i/M={Mi.min()/M:.3f} mean d_e={de.mean():.2f} frac(d_e<=1)={frac1:.2f} | A''_LP/M={a/M:.4f}"
              f" | 1/u={1/u:.4f} phi/u={phi_lb/u:.4f} | 0.564/sqrt u={0.564/u**.5:.4f} | symmetric-load formula (1/u)sum p_d E|S_d|={symf:.4f} | u*A''/min M_i={u*a/Mi.min():.2f}  [{time.time()-t:.0f}s]", flush=True)
    for u in us_prod:
        k = 2 ** u; cls = rng.permutation(np.repeat(np.arange(k), n // k))
        L = np.array([[1 if (c >> i) & 1 else -1 for c in cls] for i in range(u)]); report("product", L, u)
    for u in us_circ:
        k = 2 * u; cls = rng.permutation(np.repeat(np.arange(k), n // k))
        L = np.array([[1 if (c - i) % k < u else -1 for c in cls] for i in range(u)]); report("circular", L, u)

if __name__ == '__main__':
    run(48, 4, 1, [2, 3, 4], [2, 3, 4, 6, 8, 12])
    run(48, 3, 2, [2, 3, 4], [2, 3, 4, 6, 8])
