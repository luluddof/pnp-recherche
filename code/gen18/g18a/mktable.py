#!/usr/bin/env python3
"""Collect all *_out.txt lines into one table (table.md).  Columns: formula = exact need of (D.9) with the
paper's taper at level 0 / paper's stated bound at level 0; BB / ST = uniform need (and balanced m)."""
import re, glob, ast
rows = {}
for f in sorted(glob.glob('*_out.txt')):
    if 'kcheck' in f: continue
    for line in open(f):
        m = re.match(r"a=(\d+) L=\d+ packets=(\[.*?\]\]) M=(\d+) N=(\d+) rank=(\d+)/(\d+) log2\[sat:ST\]=(-?\d+) exh=(\w+) u=(\d+): .*?=(\{.*?\}) witness g=(\{.*?\}) \| uniform need: BB=(\w+) sat\(ST\)=(\w+) ST=(\w+)(?: \| balanced.*?BB=(\w+) ST=(\w+))?", line)
        if not m: continue
        a, pk, M, N, r, n1, idx, exh, u, F, g, bb, sat, st, bbb, stb = m.groups()
        pk = ast.literal_eval(pk); F = ast.literal_eval(F)
        name = '+'.join('[' + ''.join(map(str, p)) + ']' for p in pk)
        cnt = {}
        for p in name.split('+'): cnt[p] = cnt.get(p, 0) + 1
        name = '+'.join(p + (f'^{c}' if c > 1 else '') for p, c in cnt.items())
        rows[(int(a), int(M), name, int(u))] = (f"{r}/{n1}", F[0][1], F[0][2], bb, st, bbb or '-', stb or '-', 'exh' if exh == 'True' else 'sampled')
with open('table.md', 'w') as o:
    o.write("| a | u | packets | M | rank | (D.9) exact need, level 0 | paper bound 1+8hD+8M(a+1), level 0 | BB uniform | ST uniform | BB balanced m | ST balanced m | gens |\n|---|---|---|---|---|---|---|---|---|---|---|---|\n")
    for (a, M, name, u), v in sorted(rows.items(), key=lambda kv: (kv[0][0], kv[0][3], kv[0][1], kv[0][2])):
        o.write(f"| {a} | {u} | {name} | {M} | {v[0]} | {v[1]} | {v[2]} | {v[3]} | {v[4]} | {v[5]} | {v[6]} | {v[7]} |\n")
print(len(rows), "rows")
