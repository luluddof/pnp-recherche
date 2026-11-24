#!/usr/bin/env python3
"""G18a run5: a=3, six of eight slots filled (N=90), sampled sublattice (lower bound on ST), u=3 then u=2."""
import sys
from run1 import analyse
sl = [int(c) for c in sys.argv[1]]
analyse(3, [sl], [3, 2], sample=2500, balanced=False, stable=80, do_sat=False)
