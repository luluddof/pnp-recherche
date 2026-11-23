#!/usr/bin/env python3
"""Step 3c: growth of cutting-planes-style (MIP B&B + cuts) refutations of random 3-CNF vs PHP.
Tools: CBC (cbcbox binary) default / cuts-off, HiGHS (highspy), CaDiCaL conflicts (pysat) as resolution baseline.
"""
import json, os, random, re, subprocess, sys, time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
CBC = "/root/.local/lib/python3.13/site-packages/cbcbox/cbc_dist/bin/cbc"
TL = float(os.environ.get("TL", 60))
D = float(os.environ.get("DENS", 4.26))


def rand3(n, seed, d=D):
    from pysat.solvers import Cadical153
    rng = random.Random(seed * 1000003 + n)
    while True:
        m = round(d * n)
        cls = []
        for _ in range(m):
            vs = rng.sample(range(1, n + 1), 3)
            cls.append([v if rng.random() < .5 else -v for v in vs])
        with Cadical153(bootstrap_with=cls) as s:
            if not s.solve():
                return n, cls


def php(h):
    P = h + 1
    var = lambda i, j: i * h + j + 1
    cls = [[var(i, j) for j in range(h)] for i in range(P)]
    for j in range(h):
        for i in range(P):
            for k in range(i + 1, P):
                cls.append([-var(i, j), -var(k, j)])
    return P * h, cls


def write_lp(n, cls, path):
    L = ["Minimize", " obj: 0 x1", "Subject To"]
    for t, c in enumerate(cls):
        terms, rhs = [], 1
        for l in c:
            if l > 0:
                terms.append(f"+ x{l}")
            else:
                terms.append(f"- x{-l}"); rhs -= 1
        L.append(f" c{t}: " + " ".join(terms) + f" >= {rhs}")
    L.append("Binary")
    L += [f" x{i}" for i in range(1, n + 1)]
    L.append("End")
    open(path, "w").write("\n".join(L) + "\n")


def run_cbc(lp, mode):
    args = [CBC, lp, "-sec", str(TL), "-threads", "1"]
    if mode == "nocuts":
        args += ["-cuts", "off", "-preprocess", "off", "-heuristics", "off"]
    elif mode == "gomory":  # only Gomory (CG-type) cuts, no clique/probing/etc.
        args += ["-cuts", "off", "-gomory", "on", "-preprocess", "off", "-heuristics", "off"]
    args += ["-solve", "-quit"]
    t = time.time()
    out = subprocess.run(args, capture_output=True, text=True, timeout=TL + 60).stdout
    dt = time.time() - t
    infeas = bool(re.search(r"infeasible", out, re.I)) and "Stopped on" not in out
    m = re.search(r"Enumerated nodes:\s+(\d+)", out)
    nodes = int(m.group(1)) if m else 0
    gen = sum(int(x) for x in re.findall(r"^  [A-Za-z0-9]+\s+(\d+)\s+[\d.]+\s+\d+\s+[\d.]+\s+\S", out, re.M))
    kept = sum(int(x) for x in re.findall(r"Cut generation complete - (\d+) cuts", out))
    return dict(ok=infeas, nodes=nodes, cuts=gen, kept=kept, t=round(dt, 3))


def run_highs(lp):
    import highspy
    h = highspy.Highs(); h.setOptionValue("output_flag", False)
    h.setOptionValue("time_limit", TL); h.setOptionValue("threads", 1)
    h.readModel(lp)
    t = time.time(); h.run(); dt = time.time() - t
    st = h.modelStatusToString(h.getModelStatus())
    return dict(ok="nfeasible" in st, nodes=int(h.getInfo().mip_node_count), t=round(dt, 3), st=st)


def run_cadical(n, cls):
    from pysat.solvers import Cadical195
    with Cadical195(bootstrap_with=cls) as s:
        s.conf_budget(int(os.environ.get("CONFB", 3000000)))
        r = s.solve_limited(); st = s.accum_stats()
    return dict(ok=(r is False), conflicts=st.get("conflicts", 0))


def job(spec):
    fam, size, seed = spec
    n, cls = rand3(size, seed) if fam == "rand3" else php(size)
    lp = os.path.join(HERE, "lp", f"{fam}_{size}_{seed}.lp")
    os.makedirs(os.path.dirname(lp), exist_ok=True)
    write_lp(n, cls, lp)
    rec = dict(fam=fam, size=size, seed=seed, nvars=n, m=len(cls))
    rec["cadical"] = run_cadical(n, cls)
    rec["cbc"] = run_cbc(lp, "default")
    rec["cbc_nocuts"] = run_cbc(lp, "nocuts")
    rec["cbc_gomory"] = run_cbc(lp, "gomory")
    rec["highs"] = run_highs(lp)
    print(json.dumps(rec), flush=True)
    return rec


if __name__ == "__main__":
    fam = sys.argv[1]
    sizes = [int(x) for x in sys.argv[2].split(",")]
    seeds = range(int(sys.argv[3])) if fam == "rand3" else [0]
    specs = [(fam, s, k) for s in sizes for k in seeds]
    with Pool(4) as p:
        recs = p.map(job, specs, chunksize=1)
    json.dump(recs, open(os.path.join(HERE, f"res_{fam}_{sys.argv[2]}.json"), "w"))
