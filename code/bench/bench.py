#!/usr/bin/env python3
"""Shared SAT benchmark harness: instance families, ground truth, checking,
timing/exponent fits and a counterexample minimizer.

Candidate interface
-------------------
Python:  solve(n, clauses) -> ("SAT", assignment) | ("UNSAT", None) | ("UNKNOWN", None)
         clauses: list of lists of nonzero ints (DIMACS literals), vars 1..n.
         assignment: dict {var: bool} or iterable of literals (e.g. [1,-2,3]).
         Unassigned vars are treated as unassigned: every clause must contain a
         TRUE literal, otherwise the SAT answer is WRONG.
External: any executable reading DIMACS on stdin, printing on stdout either
         "SAT" (or "s SATISFIABLE") followed by literals (optionally on "v" lines,
         optional trailing 0), or "UNSAT"/"s UNSATISFIABLE", or "UNKNOWN".

Candidate spec strings (CLI and load_candidate()):
  path/to/file.py:func        python function
  exec:/path/to/bin [args]    external executable (shlex-split)
  builtin:unit_prop_only | builtin:brute_force

Usage
-----
  python3 bench.py --candidate cand.py:solve --suite quick [--timeout 2] [--families a,b]
  python3 bench.py --candidate cand.py:solve --minimize            # first failure in quick suite
  python3 bench.py --candidate cand.py:solve --minimize --cnf f.cnf
  python3 bench.py --candidate cand.py:solve --minimize --family rand3sat --size 40 --seed 3
     minimizer flags: --preserve-truth (keep SAT/UNSAT status), --keep-width (no literal
                      deletion / variable substitution, keeps clause widths), --min-out out.cnf
  python3 bench.py --selftest
"""
import argparse, hashlib, importlib.util, itertools, json, math, os, random, shlex
import signal, subprocess, sys, time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE_PATH = os.path.join(HERE, "truth_cache.json")
BRUTE_MAX_N = 16

# --------------------------------------------------------------------------
# DIMACS helpers
# --------------------------------------------------------------------------

def to_dimacs(n, clauses, comment=None):
    out = []
    if comment:
        out += ["c " + line for line in str(comment).splitlines()]
    out.append(f"p cnf {n} {len(clauses)}")
    out += [" ".join(map(str, c)) + " 0" for c in clauses]
    return "\n".join(out) + "\n"


def parse_dimacs(text):
    n, clauses, cur = 0, [], []
    for line in text.splitlines():
        line = line.strip()
        if not line or line[0] in "c%":
            continue
        if line[0] == "p":
            n = int(line.split()[2]); continue
        for tok in line.split():
            v = int(tok)
            if v == 0:
                clauses.append(cur); cur = []
            else:
                cur.append(v)
    if cur:
        clauses.append(cur)
    n = max([n] + [abs(l) for c in clauses for l in c])
    return n, clauses


def formula_hash(n, clauses):
    return hashlib.sha1(to_dimacs(n, clauses).encode()).hexdigest()[:16]

# --------------------------------------------------------------------------
# Generators. Each returns (n, clauses, known) where known is True/False if the
# satisfiability status is proven by construction, else None.
# --------------------------------------------------------------------------

def _rand_clause(n, rng, k=3):
    vs = rng.sample(range(1, n + 1), k)
    return [v if rng.random() < .5 else -v for v in vs]


def gen_rand3sat(n, rng, d=4.26):
    m = int(round(d * n))
    return n, [_rand_clause(n, rng) for _ in range(m)], None


def gen_planted3sat(n, rng, d=4.26):
    hidden = {v: rng.random() < .5 for v in range(1, n + 1)}
    cls = []
    while len(cls) < int(round(d * n)):
        c = _rand_clause(n, rng)
        if any(hidden[abs(l)] == (l > 0) for l in c):
            cls.append(c)
    return n, cls, True


def xor_clauses(vs, parity):
    """CNF for XOR(vs) = parity (2^(k-1) clauses)."""
    out = []
    for signs in itertools.product([0, 1], repeat=len(vs)):
        if sum(signs) % 2 != parity:
            out.append([(-v if s else v) for v, s in zip(vs, signs)])
    return out


def gf2_solvable(eqs):
    """eqs: list of (list_of_vars, parity). Gaussian elimination over GF(2)."""
    pivots = {}  # pivot bit -> (row mask, parity)
    for vs, p in eqs:
        row = 0
        for v in vs:
            row ^= 1 << v
        while row:
            hb = row.bit_length() - 1
            if hb in pivots:
                r2, p2 = pivots[hb]; row ^= r2; p ^= p2
            else:
                pivots[hb] = (row, p); break
        if row == 0 and p == 1:
            return False
    return True


def _xor_system(n, m, rng, planted=None):
    eqs = []
    for _ in range(m):
        vs = rng.sample(range(1, n + 1), 3)
        p = rng.randrange(2) if planted is None else sum(planted[v] for v in vs) % 2
        eqs.append((vs, p))
    return eqs


def gen_xor3(n, rng, d=0.92, want=None):
    """Random 3-XOR with m=d*n equations. want='sat'/'unsat' rejection-samples."""
    m = max(1, int(round(d * n)))
    for _ in range(2000):
        eqs = _xor_system(n, m, rng)
        s = gf2_solvable(eqs)
        if want is None or s == (want == "sat"):
            break
        if want == "unsat" and m < n + 4:
            m += 1  # make UNSAT reachable on tiny n
    cls = [c for vs, p in eqs for c in xor_clauses(vs, p)]
    return n, cls, s


def gen_planted_xor3(n, rng, d=1.0):
    hidden = {v: rng.randrange(2) for v in range(1, n + 1)}
    eqs = _xor_system(n, max(1, int(round(d * n))), rng, planted=hidden)
    return n, [c for vs, p in eqs for c in xor_clauses(vs, p)], True


def random_cubic_graph(nv, rng):
    while True:
        pts = [v for v in range(nv) for _ in range(3)]
        rng.shuffle(pts)
        edges = [(pts[2 * i], pts[2 * i + 1]) for i in range(len(pts) // 2)]
        if all(a != b for a, b in edges) and len({tuple(sorted(e)) for e in edges}) == len(edges):
            return edges


def gen_tseitin(nv, rng, odd=True):
    """Tseitin on a random simple 3-regular graph with nv (even) vertices; n=1.5*nv.
    Odd total charge => UNSAT, even => SAT."""
    nv += nv % 2
    edges = random_cubic_graph(nv, rng)
    inc = [[] for _ in range(nv)]
    for i, (a, b) in enumerate(edges):
        inc[a].append(i + 1); inc[b].append(i + 1)
    charge = [rng.randrange(2) for _ in range(nv)]
    if (sum(charge) % 2 == 1) != odd:
        charge[0] ^= 1
    cls = [c for v in range(nv) for c in xor_clauses(inc[v], charge[v])]
    return len(edges), cls, not odd


def gen_php(h, rng=None):
    """PHP(h+1, h): h+1 pigeons, h holes. var p(i,j)=i*h+j+1. Always UNSAT."""
    P = h + 1
    var = lambda i, j: i * h + j + 1
    cls = [[var(i, j) for j in range(h)] for i in range(P)]
    for j in range(h):
        for i in range(P):
            for k in range(i + 1, P):
                cls.append([-var(i, j), -var(k, j)])
    return P * h, cls, False


SMALL_DENSITIES = [1, 2, 3, 3.5, 4, 4.26, 4.5, 5, 6, 8, 10]


def gen_small(n, rng, seed=0):
    """'All small formulas': random 3-CNF, n in 4..10, density cycles with seed."""
    d = SMALL_DENSITIES[seed % len(SMALL_DENSITIES)]
    return n, [_rand_clause(n, rng) for _ in range(max(1, int(round(d * n))))], None


# family name -> (generator(size, rng, seed) , description of size param)
FAMILIES = {
    "small":        (lambda s, r, k: gen_small(s, r, k), "n (density cycles over seeds)"),
    "rand3sat":     (lambda s, r, k: gen_rand3sat(s, r), "n, density 4.26"),
    "planted3sat":  (lambda s, r, k: gen_planted3sat(s, r), "n, density 4.26, planted"),
    "tseitin":      (lambda s, r, k: gen_tseitin(s, r, odd=True), "graph vertices (n=1.5*V), UNSAT"),
    "tseitin_even": (lambda s, r, k: gen_tseitin(s, r, odd=False), "graph vertices (n=1.5*V), SAT"),
    "php":          (lambda s, r, k: gen_php(s), "holes h, PHP(h+1,h), UNSAT"),
    "xor3_sat":     (lambda s, r, k: gen_xor3(s, r, want="sat"), "n, m~0.92n, SAT"),
    "xor3_unsat":   (lambda s, r, k: gen_xor3(s, r, want="unsat"), "n, m~0.92n, UNSAT"),
    "planted_xor3": (lambda s, r, k: gen_planted_xor3(s, r), "n, m=n, planted"),
}


def make_instance(family, size, seed):
    rng = random.Random(f"{family}|{size}|{seed}")
    n, cls, known = FAMILIES[family][0](size, rng, seed)
    return n, cls, known


def instance_key(family, size, seed):
    return f"{family}:{size}:{seed}"

# --------------------------------------------------------------------------
# Ground truth
# --------------------------------------------------------------------------

def brute_force_sat(n, clauses):
    """Vectorised exhaustive check (n <= ~22). Returns (sat, model_dict or None)."""
    if any(len(c) == 0 for c in clauses):
        return False, None
    if n == 0:
        return True, {}
    A = np.arange(1 << n, dtype=np.uint32)
    alive = np.ones(1 << n, dtype=bool)
    for c in clauses:
        s = np.zeros(1 << n, dtype=bool)
        for l in c:
            bit = ((A >> (abs(l) - 1)) & 1).astype(bool)
            s |= bit if l > 0 else ~bit
        alive &= s
    idx = np.flatnonzero(alive)
    if len(idx) == 0:
        return False, None
    a = int(idx[0])
    return True, {v: bool((a >> (v - 1)) & 1) for v in range(1, n + 1)}


def cadical_sat(n, clauses):
    from pysat.solvers import Cadical153
    if any(len(c) == 0 for c in clauses):
        return False
    with Cadical153(bootstrap_with=clauses) as s:
        return s.solve()


def truth_of(n, clauses, known=None):
    """Returns (is_sat, source)."""
    if n <= BRUTE_MAX_N:
        s = brute_force_sat(n, clauses)[0]
        if known is not None and known != s:
            raise RuntimeError("ground-truth disagreement brute vs construction")
        return s, "brute"
    if known is not None:
        return known, "construction"
    return cadical_sat(n, clauses), "cadical"


class TruthCache:
    def __init__(self, path=CACHE_PATH):
        self.path, self.dirty = path, False
        try:
            with open(path) as f:
                self.d = json.load(f)
        except Exception:
            self.d = {}

    def get(self, key, n, clauses, known=None):
        h = formula_hash(n, clauses)
        e = self.d.get(key)
        if e is None or e.get("hash") != h:
            s, src = truth_of(n, clauses, known)
            e = self.d[key] = {"sat": s, "src": src, "hash": h, "n": n, "m": len(clauses)}
            self.dirty = True
        return e["sat"]

    def save(self):
        if self.dirty:
            tmp = self.path + ".tmp"
            with open(tmp, "w") as f:
                json.dump(self.d, f)
            os.replace(tmp, self.path)
            self.dirty = False

# --------------------------------------------------------------------------
# Candidates
# --------------------------------------------------------------------------

class CandTimeout(Exception):
    pass


def _alarm(signum, frame):
    raise CandTimeout()


def run_candidate(cand, n, clauses, timeout):
    """Returns (status, assignment, elapsed, err). status in SAT/UNSAT/UNKNOWN/TIMEOUT/ERROR."""
    t0 = time.perf_counter()
    if isinstance(cand, ExternalCandidate):
        return cand.run(n, clauses, timeout)
    old = signal.signal(signal.SIGALRM, _alarm)
    signal.setitimer(signal.ITIMER_REAL, timeout)
    try:
        res = cand(n, [list(c) for c in clauses])
        signal.setitimer(signal.ITIMER_REAL, 0)
        el = time.perf_counter() - t0
    except CandTimeout:
        return "TIMEOUT", None, time.perf_counter() - t0, None
    except Exception as e:
        signal.setitimer(signal.ITIMER_REAL, 0)
        return "ERROR", None, time.perf_counter() - t0, f"{type(e).__name__}: {e}"
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old)
    try:
        st, a = (res, None) if isinstance(res, str) else (res[0], res[1])
        st = str(st).upper()
        if st not in ("SAT", "UNSAT", "UNKNOWN"):
            return "ERROR", None, el, f"bad status {st!r}"
        return st, _norm_assign(a) if st == "SAT" else None, el, None
    except Exception as e:
        return "ERROR", None, el, f"malformed result: {e}"


def _norm_assign(a):
    if a is None:
        return {}
    if isinstance(a, dict):
        return {int(k): bool(v) for k, v in a.items()}
    return {abs(int(l)): int(l) > 0 for l in a if int(l) != 0}


class ExternalCandidate:
    def __init__(self, cmd):
        self.cmd = shlex.split(cmd)
        self.__name__ = os.path.basename(self.cmd[0])

    def run(self, n, clauses, timeout):
        t0 = time.perf_counter()
        try:
            p = subprocess.run(self.cmd, input=to_dimacs(n, clauses), capture_output=True,
                               text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            return "TIMEOUT", None, time.perf_counter() - t0, None
        el = time.perf_counter() - t0
        toks = []
        for line in p.stdout.splitlines():
            line = line.strip()
            if not line or line.startswith("c"):
                continue
            if line.startswith("s "):
                line = line[2:]
            if line.startswith("v "):
                line = line[2:]
            toks += line.split()
        if not toks:
            return "ERROR", None, el, f"no output (rc={p.returncode}) {p.stderr[-200:]}"
        st = toks[0].upper()
        st = {"SATISFIABLE": "SAT", "UNSATISFIABLE": "UNSAT"}.get(st, st)
        if st not in ("SAT", "UNSAT", "UNKNOWN"):
            return "ERROR", None, el, f"bad status {toks[0]!r}"
        if st != "SAT":
            return st, None, el, None
        try:
            return "SAT", _norm_assign([int(t) for t in toks[1:]]), el, None
        except ValueError as e:
            return "ERROR", None, el, f"bad assignment: {e}"


def load_candidate(spec):
    if spec.startswith("builtin:"):
        return BUILTINS[spec.split(":", 1)[1]]
    if spec.startswith("exec:"):
        return ExternalCandidate(spec[5:])
    if ":" in spec and spec.rsplit(":", 1)[0].endswith(".py"):
        path, fn = spec.rsplit(":", 1)
        mspec = importlib.util.spec_from_file_location(
            "cand_" + hashlib.md5(path.encode()).hexdigest()[:8], os.path.abspath(path))
        mod = importlib.util.module_from_spec(mspec)
        sys.path.insert(0, os.path.dirname(os.path.abspath(path)))
        mspec.loader.exec_module(mod)
        return getattr(mod, fn)
    if os.path.isfile(spec) and os.access(spec, os.X_OK):
        return ExternalCandidate(spec)
    raise ValueError(f"cannot load candidate {spec!r}")

# --------------------------------------------------------------------------
# Checking
# --------------------------------------------------------------------------

def satisfies(clauses, a):
    return all(any(a.get(abs(l)) is not None and a[abs(l)] == (l > 0) for l in c) for c in clauses)


def classify(status, assign, n, clauses, truth_fn):
    """Outcome: SAT_OK, UNSAT_OK, UNKNOWN, TIMEOUT, ERROR, WRONG.
    truth_fn() -> bool, only called when needed. Returns (outcome, detail)."""
    if status == "SAT":
        if satisfies(clauses, assign):
            return "SAT_OK", None
        return "WRONG", "SAT with assignment that falsifies a clause"
    if status == "UNSAT":
        if truth_fn():
            return "WRONG", "UNSAT on a satisfiable formula"
        return "UNSAT_OK", None
    return status, None


def evaluate(cand, n, clauses, timeout, truth_fn):
    st, a, el, err = run_candidate(cand, n, clauses, timeout)
    out, det = classify(st, a, n, clauses, truth_fn)
    return out, (det or err), el

# --------------------------------------------------------------------------
# Suites and reporting
# --------------------------------------------------------------------------

SUITES = {
    # family: (sizes, seeds per size)
    "quick": {
        "small":        (list(range(4, 11)), 22),
        "rand3sat":     ([10, 20, 30, 40, 60, 80, 120], 5),
        "planted3sat":  ([10, 20, 30, 40, 60, 80, 120], 5),
        "tseitin":      ([4, 6, 8, 10, 14, 20, 30], 3),
        "tseitin_even": ([4, 6, 8, 10, 14, 20, 30], 3),
        "php":          ([2, 3, 4, 5, 6, 7], 1),
        "xor3_sat":     ([8, 12, 16, 24, 32, 48], 3),
        "xor3_unsat":   ([8, 12, 16, 24, 32, 48], 3),
        "planted_xor3": ([8, 12, 16, 24, 32, 48], 3),
    },
    "full": {
        "small":        (list(range(4, 11)), 110),
        "rand3sat":     ([10, 20, 30, 40, 60, 80, 120, 160, 200, 300], 10),
        "planted3sat":  ([10, 20, 30, 40, 60, 80, 120, 160, 200, 300, 400], 10),
        "tseitin":      ([4, 6, 8, 10, 14, 20, 30, 40, 60, 80], 5),
        "tseitin_even": ([4, 6, 8, 10, 14, 20, 30, 40, 60, 80], 5),
        "php":          ([2, 3, 4, 5, 6, 7, 8, 9], 1),
        "xor3_sat":     ([8, 12, 16, 24, 32, 48, 64, 96, 128], 6),
        "xor3_unsat":   ([8, 12, 16, 24, 32, 48, 64, 96, 128], 6),
        "planted_xor3": ([8, 12, 16, 24, 32, 48, 64, 96, 128], 6),
    },
}

OUTCOMES = ["SAT_OK", "UNSAT_OK", "UNKNOWN", "TIMEOUT", "WRONG", "ERROR"]


def iter_suite(suite, families=None):
    for fam, (sizes, seeds) in SUITES[suite].items():
        if families and fam not in families:
            continue
        for size in sizes:
            for seed in range(seeds):
                yield fam, size, seed


def run_suite(cand, suite="quick", timeout=2.0, families=None, fail_dir=None,
              cache=None, verbose=True, skip_after_timeouts=2):
    """Runs the suite. Returns list of per-instance records.
    If the first `skip_after_timeouts` instances of a size all time out, the rest of
    that size and all larger sizes of the family are skipped (recorded as SKIPPED)."""
    cache = cache or TruthCache()
    recs, dead, tcount = [], set(), {}
    for fam, size, seed in iter_suite(suite, families):
        key = instance_key(fam, size, seed)
        if fam in dead:
            recs.append(dict(family=fam, size=size, seed=seed, outcome="SKIPPED")); continue
        n, cls, known = make_instance(fam, size, seed)
        truth = cache.get(key, n, cls, known)
        out, det, el = evaluate(cand, n, cls, timeout, lambda: truth)
        r = dict(family=fam, size=size, seed=seed, n=n, m=len(cls), truth=truth,
                 outcome=out, time=el, detail=det)
        recs.append(r)
        if out == "TIMEOUT":
            tcount[(fam, size)] = tcount.get((fam, size), 0) + 1
            if seed + 1 == tcount[(fam, size)] >= min(skip_after_timeouts, SUITES[suite][fam][1]):
                dead.add(fam)
        if out in ("WRONG", "ERROR", "UNKNOWN") and fail_dir:
            os.makedirs(fail_dir, exist_ok=True)
            with open(os.path.join(fail_dir, f"{out}_{fam}_{size}_{seed}.cnf"), "w") as f:
                f.write(to_dimacs(n, cls, f"{key} truth={'SAT' if truth else 'UNSAT'} "
                                          f"outcome={out} {det or ''}"))
        if verbose and out in ("WRONG", "ERROR"):
            print(f"  !! {out} on {key} (n={n}): {det}", file=sys.stderr)
    cache.save()
    return recs


def _fit(xs, ys):
    """Least squares slope, intercept, R^2."""
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    if len(xs) < 2 or np.ptp(xs) == 0:
        return float("nan"), float("nan"), float("nan")
    b, a = np.polyfit(xs, ys, 1)
    res = ys - (a + b * xs)
    tot = ((ys - ys.mean()) ** 2).sum()
    return b, a, (1 - (res ** 2).sum() / tot) if tot > 0 else 1.0


def summarize(recs, timeout=None):
    """Returns (rows, fam_summary). rows: per family+size counts; fam_summary: fits."""
    groups = {}
    for r in recs:
        groups.setdefault((r["family"], r["size"]), []).append(r)
    rows = []
    for (fam, size), rs in groups.items():
        ran = [r for r in rs if r["outcome"] != "SKIPPED"]
        cnt = {o: sum(r["outcome"] == o for r in ran) for o in OUTCOMES}
        rows.append(dict(family=fam, size=size,
                         n=round(float(np.mean([r["n"] for r in ran])), 1) if ran else None,
                         count=len(ran), skipped=len(rs) - len(ran), **cnt,
                         mean_time=float(np.mean([r["time"] for r in ran])) if ran else None))
    fams = {}
    for fam in dict.fromkeys(r["family"] for r in rows):
        fr = [r for r in rows if r["family"] == fam and r["count"] > 0]
        xs = [r["n"] for r in fr]
        ts = [max(r["mean_time"], 1e-6) for r in fr]
        slope, _, r2p = _fit(np.log(xs), np.log(ts))
        eb, _, r2e = _fit(xs, np.log(ts))
        local = [math.log(ts[i + 1] / ts[i]) / math.log(xs[i + 1] / xs[i])
                 for i in range(len(xs) - 1) if xs[i + 1] > xs[i]]
        tos = sum(r["TIMEOUT"] for r in fr) + sum(r["skipped"] for r in rows if r["family"] == fam)
        flags = []
        if tos:
            flags.append("TIMEOUTS")
        if max(ts) > 0.02 and local and local[-1] > max(4.0, 2 * max(1.0, np.median(local))):
            flags.append("SUPERPOLY?")
        if max(ts) > 0.02 and r2e > r2p + 0.03 and slope > 2:
            flags.append("EXP-FIT-BETTER")
        if any(r["WRONG"] or r["ERROR"] for r in fr):
            flags.append("WRONG!")
        fams[fam] = dict(slope=slope, r2_poly=r2p, exp_rate=eb, r2_exp=r2e,
                         last_local_slope=local[-1] if local else float("nan"),
                         flags=flags)
    return rows, fams


def print_report(rows, fams, out=sys.stdout):
    hdr = f"{'family':13} {'size':>5} {'n':>6} {'cnt':>4} {'satOK':>5} {'unsOK':>5} {'unk':>4} " \
          f"{'tmo':>4} {'WRONG':>5} {'err':>4} {'skip':>4} {'mean_t(s)':>10}"
    print(hdr, file=out); print("-" * len(hdr), file=out)
    for r in rows:
        mt = f"{r['mean_time']:.2e}" if r["mean_time"] is not None else "-"
        print(f"{r['family']:13} {r['size']:>5} {str(r['n']):>6} {r['count']:>4} {r['SAT_OK']:>5} "
              f"{r['UNSAT_OK']:>5} {r['UNKNOWN']:>4} {r['TIMEOUT']:>4} {r['WRONG']:>5} "
              f"{r['ERROR']:>4} {r['skipped']:>4} {mt:>10}", file=out)
    print(file=out)
    print(f"{'family':13} {'loglog_slope':>12} {'R2':>5} {'last_local':>10} {'exp_rate/var':>12} "
          f"{'R2exp':>6}  flags", file=out)
    for fam, s in fams.items():
        print(f"{fam:13} {s['slope']:12.2f} {s['r2_poly']:5.2f} {s['last_local_slope']:10.2f} "
              f"{s['exp_rate']:12.3f} {s['r2_exp']:6.2f}  {' '.join(s['flags'])}", file=out)
    tot = {o: sum(r[o] for r in rows) for o in OUTCOMES}
    print("\nTOTAL: " + "  ".join(f"{k}={v}" for k, v in tot.items()), file=out)

# --------------------------------------------------------------------------
# Counterexample minimizer
# --------------------------------------------------------------------------

def _simplify_assign(clauses, v, val):
    """Substitute var v := val. Returns new clause list (may contain empty clause)."""
    lit = v if val else -v
    out = []
    for c in clauses:
        if lit in c:
            continue
        out.append([l for l in c if l != -lit])
    return out


def _merge_vars(clauses, v, lit):
    """Substitute var v := lit (lit = +-u). Drops tautologies and duplicate literals."""
    out = []
    for c in clauses:
        nc = []
        for l in c:
            x = (lit if l > 0 else -lit) if abs(l) == v else l
            if x not in nc:
                nc.append(x)
        if not any(-x in nc for x in nc):
            out.append(nc)
    return out


def _compact(n, clauses):
    """Rename variables to 1..k in order of first occurrence."""
    mp = {}
    for c in clauses:
        for l in c:
            if abs(l) not in mp:
                mp[abs(l)] = len(mp) + 1
    return len(mp), [[(mp[abs(l)] if l > 0 else -mp[abs(l)]) for l in c] for c in clauses]


def make_failure_predicate(cand, timeout, targets, preserve_truth=None, allow_empty=False):
    """Returns pred(n, clauses)->bool, true iff candidate's outcome is in `targets`
    (and, with preserve_truth=True/False, the formula's SAT status equals it)."""
    stats = {"calls": 0}

    def pred(n, clauses):
        if not allow_empty and any(len(c) == 0 for c in clauses):
            return False
        stats["calls"] += 1
        tcache = {}

        def tf():
            if "t" not in tcache:
                tcache["t"] = truth_of(n, clauses)[0]
            return tcache["t"]
        if preserve_truth is not None and tf() != preserve_truth:
            return False
        out, _, _ = evaluate(cand, n, clauses, timeout, tf)
        return out in targets
    pred.stats = stats
    return pred


def ddmin(items, test):
    """Zeller's delta debugging: returns a 1-minimal sublist with test(sublist) True."""
    n = 2
    while len(items) >= 2:
        chunk = math.ceil(len(items) / n)
        subsets = [items[i:i + chunk] for i in range(0, len(items), chunk)]
        reduced = False
        for i, s in enumerate(subsets):
            if test(s):
                items, n, reduced = s, 2, True; break
        if not reduced:
            for i in range(len(subsets)):
                comp = [x for j, s in enumerate(subsets) if j != i for x in s]
                if test(comp):
                    items, n, reduced = comp, max(n - 1, 2), True; break
        if not reduced:
            if n >= len(items):
                break
            n = min(len(items), 2 * n)
    if len(items) == 1 and test([]):
        return []
    return items


def minimize(cand, n, clauses, timeout=2.0, targets=None, preserve_truth=False,
             keep_width=False, verbose=False, allow_empty=False):
    """Shrink (n, clauses) while the candidate keeps failing.
    targets: set of outcomes counted as failure (default: derived from the original
    outcome: WRONG->{WRONG}; UNKNOWN->{UNKNOWN,WRONG}; TIMEOUT->{TIMEOUT}; ERROR->{ERROR}).
    preserve_truth: only accept reductions keeping the original SAT/UNSAT status.
    keep_width: only delete clauses / rename vars (never shortens a clause).
    allow_empty: allow the empty clause in reduced formulas (off by default).
    Returns dict(n, clauses, outcome, truth, calls)."""
    t0_truth = truth_of(n, clauses)[0]
    orig, _, _ = evaluate(cand, n, clauses, timeout, lambda: t0_truth)
    if targets is None:
        targets = {"WRONG": {"WRONG"}, "UNKNOWN": {"UNKNOWN", "WRONG"},
                   "TIMEOUT": {"TIMEOUT"}, "ERROR": {"ERROR"}}.get(orig)
        if targets is None:
            raise ValueError(f"candidate does not fail on this formula (outcome {orig})")
    pred = make_failure_predicate(cand, timeout, targets,
                                  t0_truth if preserve_truth else None, allow_empty)
    if not pred(n, clauses):
        raise ValueError("failure not reproducible (nondeterministic candidate?)")
    cur = [list(c) for c in clauses]
    log = (lambda *a: print(*a, file=sys.stderr)) if verbose else (lambda *a: None)
    changed = True
    while changed:
        changed = False
        before = (len(cur), sum(map(len, cur)))
        # 1. delta-debug clause deletion
        cur = ddmin(cur, lambda cs: pred(n, cs))
        log(f"  ddmin -> {len(cur)} clauses")
        # 2. one-at-a-time clause deletion (cheap after ddmin)
        i = 0
        while i < len(cur):
            t = cur[:i] + cur[i + 1:]
            if pred(n, t):
                cur = t
            else:
                i += 1
        if not keep_width:
            # 3. variable elimination by substitution v:=True / v:=False
            for v in sorted({abs(l) for c in cur for l in c}):
                for val in (False, True):
                    t = _simplify_assign(cur, v, val)
                    if pred(n, t):
                        cur = t; break
            # 4. literal deletion
            for ci in range(len(cur)):
                j = 0
                while ci < len(cur) and j < len(cur[ci]):
                    t = [list(c) for c in cur]
                    del t[ci][j]
                    if pred(n, t):
                        cur = t
                    else:
                        j += 1
        # 5. variable merging / renaming: v := +-u for u < v
        vs = sorted({abs(l) for c in cur for l in c})
        if len(vs) <= 60:
            done = False
            for v in reversed(vs):
                for u in vs:
                    if u >= v or done:
                        break
                    for lit in (u, -u):
                        t = _merge_vars(cur, v, lit)
                        if keep_width and (len(t) != len(cur) or
                                           any(len(a) != len(b) for a, b in zip(t, cur))):
                            continue
                        if pred(n, t):
                            cur, done = t, True; break
                if done:
                    break  # restart the outer loop with the smaller formula
        # 6. drop duplicate clauses
        seen, ded = set(), []
        for c in cur:
            k = tuple(sorted(c))
            if k not in seen:
                seen.add(k); ded.append(c)
        if len(ded) < len(cur) and pred(n, ded):
            cur = ded
        if (len(cur), sum(map(len, cur))) < before:
            changed = True
    # 7. rename variables to 1..k (accept only if failure persists)
    k, comp = _compact(n, cur)
    mx = max([0] + [abs(l) for c in cur for l in c])
    if pred(k, comp):
        n, cur = k, comp
    elif pred(mx, cur):
        n = mx
    truth = truth_of(n, cur)[0]
    out, det, _ = evaluate(cand, n, cur, timeout, lambda: truth)
    return dict(n=n, clauses=cur, outcome=out, detail=det, truth=truth,
                calls=pred.stats["calls"], original_outcome=orig)


def find_first_failure(cand, suite="quick", timeout=2.0, families=None, cache=None,
                       targets=("WRONG", "ERROR", "UNKNOWN", "TIMEOUT"), truth=None):
    """Scan suite instances (smallest first within each family; WRONG preferred
    over UNKNOWN) and return (key, n, clauses, outcome) of the smallest failure."""
    cache = cache or TruthCache()
    best = None
    rank = {"WRONG": 0, "ERROR": 1, "UNKNOWN": 2, "TIMEOUT": 3}
    for fam, size, seed in iter_suite(suite, families):
        n, cls, known = make_instance(fam, size, seed)
        if best and (rank[best[3]] == 0 and len(cls) >= len(best[2])):
            continue
        tr = cache.get(instance_key(fam, size, seed), n, cls, known)
        if truth is not None and tr != truth:
            continue
        out, _, _ = evaluate(cand, n, cls, timeout, lambda: tr)
        if out in targets:
            score = (rank[out], len(cls))
            if best is None or score < (rank[best[3]], len(best[2])):
                best = (instance_key(fam, size, seed), n, cls, out)
    cache.save()
    return best

# --------------------------------------------------------------------------
# Reference candidates
# --------------------------------------------------------------------------

def unit_prop_only(n, clauses):
    """Polynomial, incomplete: unit propagation from the empty assignment.
    Conflict -> UNSAT; all clauses satisfied -> SAT; otherwise UNKNOWN."""
    a = {}
    changed = True
    while changed:
        changed = False
        for c in clauses:
            unk, nunk = None, 0
            for l in c:
                v = abs(l)
                if v in a:
                    if a[v] == (l > 0):
                        break
                else:
                    nunk += 1; unk = l
            else:
                if nunk == 0:
                    return ("UNSAT", None)
                if nunk == 1:
                    a[abs(unk)] = unk > 0; changed = True
    if all(any(abs(l) in a and a[abs(l)] == (l > 0) for l in c) for c in clauses):
        return ("SAT", {v: a.get(v, False) for v in range(1, n + 1)})
    return ("UNKNOWN", None)


def brute_force(n, clauses):
    """Exponential, complete: enumerate all 2^n assignments (pure Python)."""
    for bits in itertools.product((False, True), repeat=n):
        if all(any(bits[abs(l) - 1] == (l > 0) for l in c) for c in clauses):
            return ("SAT", {v + 1: bits[v] for v in range(n)})
    return ("UNSAT", None)


BUILTINS = {"unit_prop_only": unit_prop_only, "brute_force": brute_force}

# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def cmd_suite(cand, args):
    t0 = time.time()
    fams = args.families.split(",") if args.families else None
    fail_dir = args.fail_dir or os.path.join(HERE, "failures", _cand_name(args.candidate))
    recs = run_suite(cand, args.suite, args.timeout, fams, fail_dir=fail_dir)
    rows, fs = summarize(recs, args.timeout)
    print(f"candidate={args.candidate} suite={args.suite} timeout={args.timeout}s")
    print_report(rows, fs)
    print(f"wall time {time.time() - t0:.1f}s; failing instances saved in {fail_dir}")
    if args.out:
        with open(args.out, "w") as f:
            json.dump(dict(records=recs, rows=rows, families=fs), f, indent=1, default=str)


def cmd_minimize(cand, args):
    if args.cnf:
        with open(args.cnf) as f:
            n, cls = parse_dimacs(f.read())
        src = args.cnf
    elif args.family:
        n, cls, _ = make_instance(args.family, args.size, args.seed)
        src = instance_key(args.family, args.size, args.seed)
    else:
        fams = args.families.split(",") if args.families else None
        want = {"sat": True, "unsat": False}.get(getattr(args, "truth", None))
        hit = find_first_failure(cand, args.suite, args.timeout, fams, truth=want)
        if hit is None:
            print("no failure found in suite"); return
        src, n, cls, _ = hit
    t0 = time.time()
    targets = set(args.targets.split(",")) if args.targets else None
    res = minimize(cand, n, cls, args.timeout, targets, args.preserve_truth, args.keep_width,
                   verbose=args.verbose)
    txt = to_dimacs(res["n"], res["clauses"],
                    f"minimized from {src} (n={n}, m={len(cls)}); outcome={res['outcome']} "
                    f"({res['detail'] or ''}); truth={'SAT' if res['truth'] else 'UNSAT'}; "
                    f"preserve_truth={args.preserve_truth} keep_width={args.keep_width}; "
                    f"{res['calls']} candidate calls, {time.time() - t0:.1f}s")
    print(txt, end="")
    if args.min_out:
        with open(args.min_out, "w") as f:
            f.write(txt)
    return res


def _cand_name(spec):
    return "".join(ch if ch.isalnum() else "_" for ch in spec)[-60:]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidate", help="file.py:func | exec:CMD | builtin:NAME")
    ap.add_argument("--suite", default="quick", choices=list(SUITES))
    ap.add_argument("--families", help="comma-separated subset of: " + ",".join(FAMILIES))
    ap.add_argument("--timeout", type=float, default=2.0, help="per-instance seconds")
    ap.add_argument("--out", help="write JSON results here")
    ap.add_argument("--fail-dir", help="directory for failing instances (.cnf)")
    ap.add_argument("--minimize", action="store_true")
    ap.add_argument("--cnf", help="minimize this DIMACS file")
    ap.add_argument("--family"); ap.add_argument("--size", type=int); ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--truth", choices=["sat", "unsat"],
                    help="--minimize without --cnf/--family: only start from instances with this status")
    ap.add_argument("--targets", help="outcomes counting as failure, e.g. WRONG,UNKNOWN")
    ap.add_argument("--preserve-truth", action="store_true")
    ap.add_argument("--keep-width", action="store_true")
    ap.add_argument("--min-out", help="write minimized DIMACS here")
    ap.add_argument("--list-families", action="store_true")
    ap.add_argument("--dump", action="store_true", help="print DIMACS of --family/--size/--seed")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)
    if args.list_families:
        for k, (_, d) in FAMILIES.items():
            print(f"{k:13} {d}")
        return
    if args.dump:
        n, cls, _ = make_instance(args.family, args.size, args.seed)
        print(to_dimacs(n, cls, instance_key(args.family, args.size, args.seed)), end=""); return
    if args.selftest:
        return selftest()
    if not args.candidate:
        ap.error("--candidate required")
    cand = load_candidate(args.candidate)
    if args.minimize:
        cmd_minimize(cand, args)
    else:
        cmd_suite(cand, args)


def selftest():
    # ground truth consistency checks
    for fam in FAMILIES:
        for size in {"php": [2, 3], "tseitin": [4, 6, 8], "tseitin_even": [4, 6, 8]}.get(fam, [8, 12, 16]):
            for seed in range(3):
                n, cls, known = make_instance(fam, size, seed)
                if n <= BRUTE_MAX_N:
                    b = brute_force_sat(n, cls)[0]
                    assert b == cadical_sat(n, cls), (fam, size, seed)
                    assert known is None or known == b, (fam, size, seed)
    print("ground-truth cross-check (brute vs CaDiCaL vs construction): OK\n")
    for name in ("unit_prop_only", "brute_force"):
        a = argparse.Namespace(candidate=f"builtin:{name}", suite="quick", timeout=2.0,
                               families=None, out=None, fail_dir=None)
        cmd_suite(BUILTINS[name], a)
        print()
    for pt, kw, fams, tr in ((False, False, None, None), (True, False, "tseitin", None),
                             (True, False, "small", "unsat"), (True, True, "small", "unsat")):
        a = argparse.Namespace(cnf=None, family=None, suite="quick", families=fams, truth=tr,
                               timeout=2.0, targets=None, preserve_truth=pt, keep_width=kw,
                               min_out=None, verbose=False)
        print(f"# minimize unit_prop_only  preserve_truth={pt} keep_width={kw} families={fams} truth={tr}")
        cmd_minimize(unit_prop_only, a)
        print()


if __name__ == "__main__":
    main()
