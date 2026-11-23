import json, sys, re, statistics as st, math, glob, numpy as np
res = []
for f in sys.argv[1:]: res += json.load(open(f))
tab = {}
for h, f, s, lv, nd, t in res:
    mm = re.match(r'(r[\d.]+|ts)_[nv](\d+)_', f); key = (mm.group(1), int(mm.group(2)), h)
    tab.setdefault(key, []).append((lv if s == 'UNSAT' else float('inf'), t))
groups = sorted(set((k[0]) for k in tab)); hs = ['var','varla','pair','clause','linla']
for g in groups:
    ns = sorted(set(k[1] for k in tab if k[0] == g))
    print(f'== {g}: median leaves (#limit hits) / median sec')
    print('n    ' + ''.join(f'{h:>20}' for h in hs))
    for n in ns:
        row = f'{n:<5}'
        for h in hs:
            L = tab.get((g, n, h))
            if not L: row += f'{"-":>20}'; continue
            med = st.median([x[0] for x in L]); lim = sum(x[0] == float('inf') for x in L)
            row += f'{(str(int(med)) if med != float("inf") else "LIM")+("("+str(lim)+")" if lim else "")+" /"+format(st.median([x[1] for x in L]),".2f"):>20}'
        print(row)
    for h in hs:
        pts = [(n, math.log2(max(1, st.median([x[0] for x in tab[(g, n, h)]])))) for n in ns
               if (g, n, h) in tab and st.median([x[0] for x in tab[(g, n, h)]]) != float('inf')]
        if len(pts) >= 3:
            x = np.array([p[0] for p in pts]); y = np.array([p[1] for p in pts]); a, b = np.polyfit(x, y, 1)
            r2 = 1 - ((y - (a*x+b))**2).sum() / ((y - y.mean())**2).sum()
            print(f'  fit {h:7s}: leaves ~ 2^({a:.4f} n) = {2**a:.4f}^n, R2={r2:.3f}, ns={list(x)}')
