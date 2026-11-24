#!/usr/bin/env python3
"""G19a c4: exchangeable loads (r=2): fraction p_d of the wires at Hamming distance d, uniform over the sign types of
that distance.  Claim: A''/M = (1/u) sum_d p_d E|S_d|, E|S_d| = d*C(d-1,floor((d-1)/2))/2^(d-1) (S_d = +-1 walk)."""
import importlib.util, math, numpy as np, warnings
warnings.simplefilter('ignore')
spec = importlib.util.spec_from_file_location("c2", "c2_worst.py"); c2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(c2)
ES = lambda d: d * math.comb(d - 1, (d - 1) // 2) / 2 ** (d - 1)
rng = np.random.default_rng(0)
for u in (2, 3, 4, 5):
    C = c2.types(u, 'tern'); dist = np.abs(C).sum(1).astype(int)
    tests = [np.eye(u + 1)[d][1:] for d in range(1, u + 1)] + [rng.dirichlet(np.ones(u)) for _ in range(3)]
    tests.append(np.array([math.comb(u, d) / (2 ** u - 1) for d in range(1, u + 1)]))   # product-uniform loads conditioned on d>=1
    for p in tests:
        w = np.array([p[d - 1] / (dist == d).sum() for d in dist])
        v = c2.V(w, C, u); f = sum(p[d - 1] * ES(d) for d in range(1, u + 1)) / u
        print(f"u={u} p={np.round(p,3)}: LP value {v:.6f}  formula {f:.6f}  {'OK' if abs(v-f)<1e-7 else 'MISMATCH'}")
