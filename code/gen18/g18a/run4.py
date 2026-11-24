#!/usr/bin/env python3
"""G18a run4: does ONE packet violating the L2 hypothesis restore the black-box depth?"""
from run1 import analyse
analyse(2, [[0, 1], [0, 1], [0, 2]], [1, 2])
analyse(2, [[0, 1], [0, 1], [0, 1], [0, 2]], [1, 2], balanced=False)
analyse(2, [[0, 1, 2, 3], [0, 1], [0, 1]], [1, 2], balanced=False)
