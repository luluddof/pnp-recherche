"""G7a core. Commutative (Hadamard-diagonal) class, one point per direction.
P_z(t) = 1 + t R_z(t), deg R_z <= d-1  (P_z = chi(z) g_z).
Edge (z, z^k) good iff (P_z + P_{z^k})(tau_k) = 0  <=>  tau_k (R_z + R_{z^k})(tau_k) = -2.
Vertices z in {0,1}^l as ints.  Feasibility of a good set: linear system, mod prime with random tau (generic),
or exactly over Q with given tau."""
import itertools, random
import flint
P = 2**31 - 1

def edges(l):
    return [(z, z | (1 << k), k) for z in range(1 << l) for k in range(l) if not z >> k & 1]

def system(l, d, good, taus, mod=P):
    """rows for good edges; unknowns R_z coefficients (d per vertex) in monomial basis."""
    n = (1 << l) * d
    rows = []
    for (z, w, k) in good:
        r = [0] * (n + 1)
        tk = taus[k]
        pw = tk  # tau * tau^j
        for j in range(d):
            r[z * d + j] = (r[z * d + j] + pw) % mod
            r[w * d + j] = (r[w * d + j] + pw) % mod
            pw = pw * tk % mod
        r[n] = (-2) % mod
        rows.append(r)
    return rows, n

def feasible(l, d, good, taus=None, seed=1, mod=P):
    if not good: return True
    if taus is None:
        rng = random.Random(seed); taus = [rng.randrange(2, mod - 1) for _ in range(l)]
    rows, n = system(l, d, good, taus, mod)
    A = flint.nmod_mat([r[:n] for r in rows], mod)
    Ab = flint.nmod_mat(rows, mod)
    return A.rank() == Ab.rank()

def feasible_Q(l, d, good, taus):
    """exact over Q, taus rationals (fmpq)"""
    n = (1 << l) * d
    rows = []
    for (z, w, k) in good:
        r = [flint.fmpq(0)] * (n + 1)
        tk = flint.fmpq(taus[k].numerator, taus[k].denominator); pw = tk
        for j in range(d):
            r[z * d + j] += pw; r[w * d + j] += pw; pw *= tk
        r[n] = flint.fmpq(-2); rows.append(r)
    A = flint.fmpq_mat(len(rows), n, [x for r in rows for x in r[:n]])
    Ab = flint.fmpq_mat(len(rows), n + 1, [x for r in rows for x in r])
    return A.rank() == Ab.rank()

def subcubes(l, m):
    """all m-dim subcubes as lists of edge indices"""
    E = edges(l); idx = {(z, k): i for i, (z, w, k) in enumerate(E)}
    out = []
    for S in itertools.combinations(range(l), m):
        Sm = sum(1 << s for s in S)
        rest = [j for j in range(l) if j not in S]
        for bits in range(1 << len(rest)):
            base = sum(1 << rest[i] for i in range(len(rest)) if bits >> i & 1)
            es = []
            for sub in range(1 << m):
                z = base | sum(1 << S[i] for i in range(m) if sub >> i & 1)
                for k in S:
                    if not z >> k & 1: es.append(idx[(z, k)])
            out.append(es)
    return out
