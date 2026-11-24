import json, glob, math, collections
import numpy as np
rows = []
for f in glob.glob('res_d*.jsonl') + glob.glob('r2_d*.jsonl'):
    for l in open(f):
        l = l.strip()
        if l.startswith('{'): rows.append(json.loads(l))
by = collections.defaultdict(list)
for r in rows: by[r['d']].append(r)
print("d | n range | slope log2(#IS) | naive | smart(branch+comp+deg2) | frontierDP total | (ln d)^2/(2d ln2) | log2(lambda_d) | slope_smart*d/log2(d)^2")
for d in sorted(by):
    R = by[d]
    def slope(key, nmin=30):
        pts = collections.defaultdict(list)
        for r in R:
            if r.get(key) and r['n'] >= nmin:
                pts[r['n']].append(math.log2(r[key]) if key != 'log2IS' else r[key])
        ns = sorted(pts)
        if len(ns) < 3: return float('nan'), ''
        y = [np.mean(pts[n]) for n in ns]
        a, b = np.polyfit(ns, y, 1)
        return a, f"{ns[0]}-{ns[-1]}"
    sI, rg = slope('log2IS'); sN, rn = slope('naive'); sS, rs = slope('smart'); sD, rd = slope('dp_total')
    # lambda_d: root of x^{d+1}=x^d+1
    lo, hi = 1.0, 2.0
    for _ in range(100):
        mid = (lo+hi)/2
        if mid**(d+1) - mid**d - 1 > 0: hi = mid
        else: lo = mid
    pred = math.log(d)**2/(2*d*math.log(2))
    print(f"{d:2d} | {rs:7s} | {sI:.3f} | {sN:.3f} ({rn}) | {sS:.3f} | {sD:.3f} ({rd}) | {pred:.3f} | {math.log2(lo):.3f} | {sS*d/math.log2(d)**2:.2f}")
