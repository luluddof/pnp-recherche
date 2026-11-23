#!/usr/bin/env python3
# external candidate example: reads DIMACS on stdin, uses CaDiCaL
import sys
from pysat.solvers import Cadical153
sys.path.insert(0, '/tmp/claude-0/-home-user/2d12e629-e6d2-542b-bfc0-8b766f57383f/scratchpad/bench')
from bench import parse_dimacs
n, cls = parse_dimacs(sys.stdin.read())
with Cadical153(bootstrap_with=cls) as s:
    if s.solve():
        print("SAT"); print(" ".join(map(str, s.get_model() or [])), "0")
    else:
        print("UNSAT")
