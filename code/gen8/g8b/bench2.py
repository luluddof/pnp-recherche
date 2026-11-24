import sys, math, time, json
import networkx as nx
import count_is as C
from bench import limited_naive
d = int(sys.argv[1])
for n in range(40, 121, 10):
    if n*d % 2 or d >= n: continue
    for seed in range(3):
        G = nx.random_regular_graph(d, n, seed=7000*d + 10*n + seed)
        nn, adj = C.to_masks(G)
        t = time.time()
        cnt, w_sb = C.smart_branch(nn, adj)
        c1, w_nb = limited_naive(nn, adj, 4_000_000) if n <= 80 else (None, None)
        if n <= 70:
            c3, w_dp, mx = C.frontier_dp(nn, adj); assert c3 == cnt
        else:
            w_dp = mx = None
        print(json.dumps(dict(d=d, n=n, seed=seed, log2IS=math.log2(cnt), naive=w_nb, smart=w_sb, dp_total=w_dp, dp_max=mx, sec=time.time()-t)), flush=True)
        if time.time()-t > 200: sys.exit()
