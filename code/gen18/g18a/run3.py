#!/usr/bin/env python3
"""G18a run3: a=2 mixes, a=3 density census, a=3 full packet (sampled sublattice), K-doubling check."""
import sys
from run1 import analyse
which = sys.argv[1]
if which == "mix":
    for k in (2, 3): analyse(2, [[0, 1]] * k, [1, 2])
    for k in (2, 3): analyse(2, [[0, 2]] * k, [1, 2])
    analyse(2, [[0, 1], [0, 2]], [1, 2]); analyse(2, [[0, 1], [0, 3]], [1, 2]); analyse(2, [[0, 1], [2, 3]], [1, 2])
    analyse(2, [[0, 1], [0]], [1, 2]); analyse(2, [[0, 1, 2]] * 2, [1, 2])
    analyse(2, [[0, 1], [0, 2], [0, 3]], [1, 2]); analyse(2, [[0, 1, 2, 3], [0, 1]], [1, 2]); analyse(2, [[0, 1, 2, 3], [0]], [1, 2])
    analyse(2, [[0, 1, 2]] * 3, [1, 2]); analyse(2, [[0, 1, 2, 3]] * 3, [1, 2], balanced=False)
if which == "a3single":
    for sl in ([0], [0, 1], [0, 4], [0, 1, 2], [0, 1, 3], [0, 2, 5], [0, 1, 2, 3], [0, 1, 3, 6], [0, 2, 4, 6], [0, 1, 4, 5],
               [0, 1, 2, 3, 4], [0, 1, 2, 4, 7], [0, 2, 3, 5, 6]):
        analyse(3, [sl], [1, 2, 3], balanced=False)
if which == "a3multi":
    for k in (2, 3, 4): analyse(3, [[0, 1]] * k, [1, 2, 3], balanced=False)
    analyse(3, [[0, 1, 3]] * 2, [1, 2, 3], balanced=False)
    analyse(3, [[0, 2, 5], [0, 1, 3]], [1, 2, 3], balanced=False)
    analyse(3, [[0, 1, 3, 6], [0, 1]], [1, 2, 3], balanced=False)
if which == "a3full":
    analyse(3, [[0, 1, 2, 3, 4, 5, 6, 7]], [3, 2, 1], sample=3000, balanced=False, stable=100, do_sat=False)
if which == "kcheck":
    analyse(2, [[0, 1]] * 3, [1, 2]); analyse(2, [[0, 2]] * 3, [1, 2]); analyse(2, [[0, 1, 2, 3]], [1, 2]); analyse(3, [[0, 1, 3]], [1, 2, 3], balanced=False)
