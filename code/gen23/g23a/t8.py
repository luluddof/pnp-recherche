from orbit import *; from groups import catalogue
for name, G in catalogue(8).items():
    run(name, G, [('S', frozenset([0])), ('T', frozenset([0]))], caps=[4, 4], timeout=60)
