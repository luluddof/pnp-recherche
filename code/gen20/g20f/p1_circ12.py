import sys; sys.argv=['x']
exec(open('p1_circ.py').read().split('# (1d) brute force')[0])
exec('def circ_exact'+open('p1_circ.py').read().split('def circ_exact')[1].split('G=nx.random')[0])
import networkx as nx, time
for seed in (11,12):
    G=nx.random_regular_graph(4,48,seed=seed)
    t0=time.time(); M,Mi,lp,ip=circ_exact(G,12,seed)
    print(f'graph seed {seed} u=12: M/u={M/12} min M_i={Mi.min()} A_LP={lp:.3f} A_int={ip} ({time.time()-t0:.0f}s)',flush=True)
