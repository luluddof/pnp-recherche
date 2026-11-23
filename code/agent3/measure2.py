import statistics, measure as M
print("== random 3-XOR, m=1.5n (unsat) ==\nn | CDCL_conflicts DRUP_lines | GaussOps")
for n in [20, 40, 60, 80, 100]:
    res=[]; seed=0
    while len(res)<3:
        cls, rows = M.rand3xor(n, int(1.5*n), 1000*n+seed); seed+=1
        g = M.gauss_ops(rows)
        if g is None: continue
        sat, conf, plen = M.cdcl(cls); assert not sat
        res.append((conf, plen, g))
    print(n,'|',*[statistics.median(r[k] for r in res) for k in range(3)], flush=True)
print("== random 3-CNF (unsat only) ==\nn ratio | unsat/tried | med conflicts | med DRUP")
for ratio in [4.26, 6.0]:
    for n in [50, 100, 150, 200]:
        res=[]; tried=0
        while len(res)<5 and tried<20:
            cls = M.rand3cnf(n, ratio, 7*n+tried); tried+=1
            sat, conf, plen = M.cdcl(cls)
            if not sat: res.append((conf, plen))
        if res: print(n, ratio,'|','%d/%d'%(len(res),tried),'|',statistics.median(r[0] for r in res),statistics.median(r[1] for r in res), flush=True)
