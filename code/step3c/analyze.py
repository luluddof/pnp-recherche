#!/usr/bin/env python3
import json, sys, math, glob
import numpy as np

def load(path):
    return [json.loads(l) for l in open(path) if l.startswith("{")]

def r2(x, y):
    A = np.vstack([x, np.ones_like(x)]).T
    c, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ c
    ss = ((y - y.mean()) ** 2).sum()
    return c[0], 1 - ((y - yhat) ** 2).sum() / ss if ss > 0 else float("nan")

METRICS = [("cadical", "conflicts"), ("cbc", "nodes"), ("cbc", "cuts"), ("cbc", "t"),
           ("cbc_nocuts", "nodes"), ("cbc_gomory", "nodes"), ("cbc_gomory", "cuts"),
           ("highs", "nodes"), ("highs", "t")]

for path in sys.argv[1:]:
    recs = load(path)
    sizes = sorted({r["size"] for r in recs})
    print(f"\n== {path}  (median over instances; '*k' = k timeouts, counted at their censored value)")
    print("size  #  " + "  ".join(f"{a}.{b}" for a, b in METRICS))
    med = {}
    for s in sizes:
        rs = [r for r in recs if r["size"] == s]
        row = []
        for a, b in METRICS:
            v = [r[a][b] for r in rs]
            to = sum(1 for r in rs if not r[a]["ok"])
            m = float(np.median(v)); med.setdefault((a, b), []).append(m)
            row.append(f"{m:.3g}" + (f"*{to}" if to else ""))
        print(f"{s:4d} {len(rs):2d}  " + "  ".join(row))
    xs = np.array(sizes, float)
    if recs[0]["fam"] == "php":
        xs = np.array([recs[[r['size'] for r in recs].index(s)]["nvars"] for s in sizes], float)
    print("fits on log(median+1): exp = vs n(vars), poly = vs log n")
    for k, v in med.items():
        y = np.log(np.array(v) + 1)
        se, re_ = r2(xs, y); sp, rp = r2(np.log(xs), y)
        win = "exp" if re_ > rp else "poly"
        print(f"  {k[0]}.{k[1]:10s} exp: rate={se:.3f}/var R2={re_:.3f} | poly: deg={sp:.2f} R2={rp:.3f} -> {win}")
