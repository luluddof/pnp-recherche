import sys, json, numpy as np
from spread import *; from load import *
out = []
def rep(name, d): d['fam'] = name; print(f"{name:14s}", fmt(d), flush=True); out.append({k: (float(v) if not isinstance(v, str) else v) for k, v in d.items()})
for f in ['exact_4_5_3.json', 'exact_5_7_4.json', 'exact_5_6_6.json', 'exact_6_7_7.json']:
    R, C, p, q = load_g19(f); rep('g19 ' + f[6:-5], analyze(R, C))
R, C, p, q = load_g24g(); rep('g24g 6_12_3', analyze(R, C))
for f in ['shape_6_7_7_S7T7U0.json', 'shape_7_8_8_S8T8U0.json', 'shape_8_8_8_S8T8U0.json']:
    R, C, p, q = load_mono(f); rep('mono ' + f[6:12], analyze(R, C))
for n in [2, 3, 4, 5]:
    R, C = cd(n); rep(f'CD n={n}', analyze(R, C))
for r, k in [(3, 2), (4, 2), (5, 2)]:
    R, C, S = hrubes(r, k); rep(f'Hrubes r={r}k={k}', analyze(R, C))
json.dump(out, open('table1.json', 'w'), indent=0)
