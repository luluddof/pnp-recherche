#!/usr/bin/env python3
"""G19a c5: class-level LP value of the staircase (= circular labelling folded by the antipode): V = 1/u."""
import importlib.util, numpy as np, warnings
warnings.simplefilter('ignore')
spec = importlib.util.spec_from_file_location("c2", "c2_worst.py"); c2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(c2)
for u in range(2, 9):
    C = np.array([[1.0 if c >= i else -1.0 for i in range(u)] for c in range(u)])
    print(f"staircase u={u}: V = {c2.V(np.ones(u) / u, C, u):.6f}   1/u = {1/u:.6f}")
