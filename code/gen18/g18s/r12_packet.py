#!/usr/bin/env python3
"""G18S r12: full antipodal packet (one wire per slot x: row end at corner x, column end at ~x; M = 2^u), majority cut."""
import math
for u in (3, 5, 7, 9, 11):
    P = 2 ** u; full = P - 1; wires = [(x, x ^ full) for x in range(P)]
    maj = lambda x: bin(x).count("1") * 2 > u
    sel = [(p, q) if maj(p) else (q, p) for (p, q) in wires]   # (C1-end, C2-end); every wire is cut
    cs = []
    for d in range(u):
        c = 0
        for g in (0, 1):
            c += min(sum(1 for a, b in sel if (a >> d) & 1 == g), sum(1 for a, b in sel if (b >> d) & 1 == g))
        cs.append(c)
    w = len(sel); val = w - max(cs); pred = P * math.comb(u - 1, (u - 1) // 2) // 2 ** (u - 1)
    print(f"u={u}: M={P}  A' majority cut: log2 N >= {val} = {val/P:.4f} M  (formula M*C(u-1,(u-1)/2)/2^(u-1) = {pred}; *sqrt(u)/M = {val/P*math.sqrt(u):.3f})   Thm A: {2}")
