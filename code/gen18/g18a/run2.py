#!/usr/bin/env python3
"""G18a run2: multi-packet scaling at a=2 (does the sparse gain survive many packets?)."""
import sys
from run1 import analyse
which = sys.argv[1]
if which == "sp01":
    for k in range(1, 7): analyse(2, [[0, 1]] * k, [1, 2])
if which == "sp02":
    for k in range(1, 6): analyse(2, [[0, 2]] * k, [1, 2])
if which == "sp0":
    for k in range(1, 7): analyse(2, [[0]] * k, [1, 2])
if which == "full":
    for k in range(1, 4): analyse(2, [[0, 1, 2, 3]] * k, [1, 2])
if which == "mix":
    analyse(2, [[0, 1], [0, 2]], [1, 2]); analyse(2, [[0, 1], [0, 3]], [1, 2]); analyse(2, [[0, 1], [2, 3]], [1, 2])
    analyse(2, [[0, 1], [0]], [1, 2]); analyse(2, [[0, 1], [0], [0]], [1, 2])
    analyse(2, [[0, 1, 2]] * 2, [1, 2]); analyse(2, [[0, 1, 2]] * 3, [1, 2])
    analyse(2, [[0, 1], [0, 2], [0, 3]], [1, 2]); analyse(2, [[0, 1], [0, 2], [0, 3], [1, 2]], [1, 2])
