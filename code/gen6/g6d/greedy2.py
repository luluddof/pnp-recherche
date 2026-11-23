import sys, random, greedy
l, d, trials = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
data = greedy.build(l, d); V, edges, rows, n = data
rng = random.Random(7); best = None
for tr in range(trials):
    order = list(range(len(edges))); rng.shuffle(order)
    order.sort(key=lambda e: (edges[e][2] >= d))          # exact directions first
    bad = greedy.greedy(l, d, order, data)
    per = [0]*l
    for e in bad: per[edges[e][2]] += 1
    if best is None or len(bad) < best[0]: best = (len(bad), per)
print(f"l={l} d={d} K-first greedy best total={best[0]} per_dir={best[1]} conj={(l-d)*2**(l-d-1)}")
