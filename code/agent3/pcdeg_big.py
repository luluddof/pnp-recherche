import sys, time, idea
idea.P = 3
for nv in [12, 14, 16]:
    ax, ne = idea.tseitin_axioms(nv, 1)
    for D in [4, 5]:
        t = time.time()
        E = idea.pc_closure(ax, ne, D)
        print(nv, ne, 'D=%d' % D, 'refuted' if frozenset() in E.rows else 'not refuted', 'dimV=%d' % len(E.rows), '%.0fs' % (time.time() - t), flush=True)
        if frozenset() in E.rows: break
