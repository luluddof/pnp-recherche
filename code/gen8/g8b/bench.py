import sys, math, time, json
import networkx as nx
import count_is as C

class Abort(Exception):
    pass

def limited_naive(n, adj, cap):
    work = [0]
    def rec(S):
        work[0] += 1
        if work[0] > cap:
            raise Abort
        best, bd = -1, 0
        for v in C.bits(S):
            dv = bin(adj[v] & S).count("1")
            if dv > bd:
                best, bd = v, dv
        if bd == 0:
            return 1 << bin(S).count("1")
        return rec(S & ~(1 << best)) + rec(S & ~(1 << best) & ~adj[best])
    try:
        return rec((1 << n) - 1), work[0]
    except Abort:
        return None, None

if __name__ != "__main__":
    pass
if __name__ == '__main__':
    d = int(sys.argv[1]) if __name__ == "__main__" else None
    out = []
    for n in [20, 30, 40, 50, 60]:
        if n * d % 2 or d >= n:
            continue
        for seed in range(3):
            G = nx.random_regular_graph(d, n, seed=1000 * d + 10 * n + seed)
            nn, adj = C.to_masks(G)
            t = time.time()
            cnt, w_dp, mx = C.frontier_dp(nn, adj)
            c2, w_sb = C.smart_branch(nn, adj)
            assert c2 == cnt
            c1, w_nb = limited_naive(nn, adj, 2_000_000)
            if c1 is not None:
                assert c1 == cnt
            rec = dict(d=d, n=n, seed=seed, log2IS=math.log2(cnt), naive=w_nb,
                       smart=w_sb, dp_total=w_dp, dp_max=mx, sec=time.time() - t)
            out.append(rec)
            print(json.dumps(rec), flush=True)
