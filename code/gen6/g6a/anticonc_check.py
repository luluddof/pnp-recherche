# check E|X-a| >= 0.216*min(V, sqrt V) for sums of independent Bernoullis, worst integer/real a
import random, math
rng = random.Random(1); worst = 9
for trial in range(4000):
    k = rng.randint(1, 40)
    ps = [rng.choice([rng.random(), rng.random()**4, 1 - rng.random()**4, 1/3, 2/3, 0.5]) for _ in range(k)]
    dist = [1.0]
    for p in ps:
        nd = [0.0]*(len(dist)+1)
        for i, q in enumerate(dist): nd[i] += q*(1-p); nd[i+1] += q*p
        dist = nd
    V = sum(p*(1-p) for p in ps)
    if V < 1e-9: continue
    # E|X-a| is convex piecewise linear in a, minimised at a median (an integer)
    m = min(sum(q*abs(i-a) for i, q in enumerate(dist)) for a in range(len(dist)))
    worst = min(worst, m/(min(V, math.sqrt(V))))
print("min over trials of E|X-a| / min(V,sqrt V) =", round(worst, 4), "(claimed >= 0.216)")
