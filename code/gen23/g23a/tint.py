from orbit import *; from groups import catalogue
import time
G = catalogue(16)['Q8xZ2']; mul, inv = G[0], G[1]
H2 = normal_subgroups(mul, inv, 2)[0]
sp_ = [('S', frozenset([0])), ('S', H2), ('T', frozenset([0])), ('T', H2)]
t = time.time(); print(run('Q8xZ2', G, sp_, caps=[7, 3, 7, 3], timeout=10, tries=5), time.time() - t)
