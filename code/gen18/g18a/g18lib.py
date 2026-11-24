#!/usr/bin/env python3
"""G18a library: exact 2-adic lattice computations for App. D toy instances.
All arithmetic is in Z/2^K (images of 2-adic integers); valuations below K_eff are exact."""
import itertools, math, sys, time
from fractions import Fraction as Fr
from math import comb
import numpy as np

import os
K = int(os.environ.get('G18K', '2048'))
MASK = (1 << K) - 1
BIG = 10 ** 9
TH = K // 2            # entries divisible by 2^TH are treated as zero (junk from truncated quotients)
THM = (1 << TH) - 1

def v2i(x):
    return (x & -x).bit_length() - 1 if x else BIG
def v2f(x):
    x = Fr(x)
    if x == 0: return BIG
    n, d, k = abs(x.numerator), x.denominator, 0
    while n % 2 == 0: n //= 2; k += 1
    while d % 2 == 0: d //= 2; k -= 1
    return k
def tomod(fr):
    fr = Fr(fr)
    return (fr.numerator * pow(fr.denominator, -1, 1 << K)) & MASK
VAL = np.frompyfunc(v2i, 1, 1)

# ---------------- polynomials mod 2^K (lists, low degree first), Kronecker multiplication
WB = (2 * K + 64) // 8
def pack(c):
    return int.from_bytes(b''.join(x.to_bytes(WB, 'little') for x in c), 'little')
def unpack(x, n):
    b = x.to_bytes(n * WB, 'little')
    return [int.from_bytes(b[i * WB:(i + 1) * WB], 'little') & MASK for i in range(n)]
def pmul(A, B):
    return unpack(pack(A) * pack(B), len(A) + len(B) - 1)

class Setup:
    def __init__(self, a):
        self.a = a; self.L = L = 2 ** a; self.S = S = list(range(-L, L)); self.al = 2 ** (a + 1)
        cp0 = 1
        for q in S:
            if q != 0: cp0 *= -q
        self.e = v2i(abs(cp0)); assert self.e == 2 * L - 2 - a
        odd = Fr(cp0, 2 ** self.e)
        self.wh = {}                      # wh[d] = 2^e * w_d  (2-integral coefficients), mod 2^K
        for d in range(-L + 1, L):
            p = [1]
            for q in S:
                if q != d:
                    p = [(-q) * p[0]] + [p[i - 1] - q * p[i] for i in range(1, len(p))] + [p[-1]]
            self.wh[d] = [tomod(Fr(c) / odd) for c in p]

# ---------------- triangular Z_(2)-basis of a lattice given by generators (vectors mod 2^K)
def insert(basis, v):
    n = len(v); c = 0; changed = False
    while True:
        while c < n and (v[c] & THM) == 0:
            v[c] = 0; c += 1
        if c == n: return changed
        b = basis.get(c)
        if b is None:
            basis[c] = v; return True
        vv = v2i(v[c]); vb = v2i(b[c])
        if vv < vb:
            basis[c] = v; v, b = b, v; vb = vv; changed = True
        f = ((v[c] >> vb) * pow(b[c] >> vb, -1, 1 << K)) & MASK
        v = [(x - f * y) & MASK for x, y in zip(v, b)]
        assert v[c] == 0

def packet_lattice(st, slots, sample=None, rng=None, stable=400, log=None):
    """triangular basis (in coordinates c' = 2^{sigma e} * monomial coefficients) of Lam_ST of one packet.
    sample=None: exhaustive over all X and all row subsets R.  sample=n: random (X,R) until `stable`
    consecutive generators do not change the lattice (gives a SUBLATTICE = rigorous lower bound only)."""
    sig = len(slots); L = st.L; e = st.e; n = sig * (2 * L - 1) + 1
    basis = {}
    def forms(X):
        fs = []
        for x in slots:
            f = [0] * (2 * L)
            for k, y in enumerate(slots):
                if X[k]: f = [(u + w) & MASK for u, w in zip(f, st.wh[y - x])]
            fs.append(f)
        return fs
    def emit(p, r):
        v = [(c << ((sig - r) * e)) & MASK for c in p] + [0] * (n - len(p))
        return insert(basis, v)
    emit([1], 0)
    ngen = 0
    if sample is None:
        seen = set()
        for X in itertools.product((0, 1), repeat=sig):
            if not any(X): continue
            fs = forms(X)
            prods = {0: [1]}
            for R in range(1, 1 << sig):
                low = (R & -R).bit_length() - 1
                p = pmul(prods[R & (R - 1)], fs[low]); prods[R] = p
                key = (bin(R).count('1'), hash(tuple(p)))
                if key in seen: continue
                seen.add(key); emit(p, bin(R).count('1')); ngen += 1
        return basis, ngen, True
    else:
        quiet = 0
        while quiet < stable and ngen < sample:
            X = [rng.randint(0, 1) for _ in range(sig)]
            if not any(X): continue
            fs = forms(X)
            order = list(range(sig)); rng.shuffle(order)
            p = [1]; r = 0; ch = False
            for x in order:                       # a random chain of nested R's
                p = pmul(p, fs[x]); r += 1
                ch |= emit(p, r); ngen += 1
            quiet = 0 if ch else quiet + 1
            if log and ngen % 2000 < sig: print(f"      .. sampled {ngen} gens, rank {len(basis)}, quiet {quiet}", file=log, flush=True)
        return basis, ngen, False

def product_lattice(b1, n1, b2, n2):
    basis = {}; n = n1 + n2 - 1
    r2 = [pack(v[:n2]) for v in b2.values()]
    for v in b1.values():
        pv = pack(v[:n1])
        for w in r2:
            insert(basis, unpack(pv * w, n))
    return basis

# ---------------- coordinates: Phi(F) = coefficients of F(t+alpha z), all t in S, j=0..N
def taylor_mats(st, N):
    mats = {}
    for t in st.S:
        A = np.empty((N + 1, N + 1), dtype=object)
        for k in range(N + 1):
            for j in range(N + 1):
                A[k, j] = (comb(k, j) * t ** (k - j) * st.al ** j) if j <= k else 0
        mats[t] = A
    return mats
def phi(st, basis, N, scale_bits, mats):
    """rows of the basis in Phi coordinates (dict t -> r x (N+1) object array), exact mod 2^(K-scale_bits)."""
    rows = np.array([v + [0] * (N + 1 - len(v)) for v in basis.values()], dtype=object)
    out = {}
    low = (1 << scale_bits) - 1
    for t in st.S:
        P = rows.dot(mats[t])
        flat = P.ravel()
        for i in range(flat.size):
            x = flat[i] & MASK
            assert x & low == 0, "kernel polynomial not integral on a disc (Lemma D.4 violated?)"
            flat[i] = x >> scale_bits
        out[t] = flat.reshape(P.shape)
    return out

def saturate(Mx):
    """Mx: r x n object array of ints mod 2^K whose rows span a lattice inside Z_(2)^n.  Returns rows of a
    basis of the SATURATION (Q-span intersected with Z_(2)^n), the pivot valuations, and the max loss."""
    R = Mx.copy(); r, n = R.shape; used = np.zeros(r, bool); vals = []; out = []
    for _ in range(r):
        rem = np.flatnonzero(~used)
        V = VAL(R[rem]).astype(np.int64); m = int(V.min())
        if m >= TH: break
        ii, c = np.argwhere(V == m)[0]; i = rem[ii]
        inv = pow(R[i, c] >> m, -1, 1 << K)
        for k in rem:
            if k == i or R[k, c] == 0: continue
            f = ((R[k, c] >> m) * inv) & MASK
            R[k] = (R[k] - f * R[i]) & MASK
        used[i] = True; vals.append(m)
        out.append(np.array([x >> m for x in R[i]], dtype=object))
    return np.array(out, dtype=object), vals

def bb_lattice(st, N, mats):
    """saturated basis of Lam_BB in Phi coordinates (all nodes)."""
    A = np.concatenate([mats[t] for t in st.S], axis=1)
    A = np.array([[x & MASK for x in row] for row in A], dtype=object)
    B, vals = saturate(A)
    assert len(vals) == N + 1, "raise K"
    return B, max(vals)

def smith_pivots(D, tgt):
    """D: r x n data matrix (object ints mod 2^K), tgt: list of r ints.  Local Smith reduction carrying the
    target.  Returns list of (pivot valuation, valuation of transformed target) and leftover targets."""
    R = D.copy(); w = list(tgt); r, n = R.shape
    ur = np.zeros(r, bool); uc = np.zeros(n, bool); piv = []
    while True:
        rem = np.flatnonzero(~ur); remc = np.flatnonzero(~uc)
        if rem.size == 0 or remc.size == 0: break
        V = VAL(R[np.ix_(rem, remc)]).astype(np.int64); m = int(V.min())
        if m >= TH: break
        ii, cc = np.argwhere(V == m)[0]; i = rem[ii]; c = remc[cc]
        inv = pow(R[i, c] >> m, -1, 1 << K)
        for k in rem:
            if k == i or R[k, c] == 0: continue
            f = ((R[k, c] >> m) * inv) & MASK
            R[k] = (R[k] - f * R[i]) & MASK; w[k] = (w[k] - f * w[i]) & MASK
        R[i, remc] = 0                         # column operations (integral since m is the global minimum)
        ur[i] = True; uc[c] = True; piv.append((m, v2i(w[i])))
    left = [v2i(w[k]) for k in np.flatnonzero(~ur)]
    return piv, left

def need_from(piv, left, p=1):
    if any(l < p for l in left): return None       # (junk has valuation >= TH)       # target not determined by the data at all
    return max([v + p - vw for v, vw in piv if vw < p], default=0)

class Instance:
    def __init__(self, st, PhiRows, N, u, label=""):
        """PhiRows: dict t -> r x (N+1) array, or a single r x 2L(N+1) array (all nodes, order st.S)."""
        self.st = st; self.N = N; self.u = u
        if isinstance(PhiRows, dict):
            self.full = np.concatenate([PhiRows[t] for t in st.S], axis=1)
        else: self.full = PhiRows
        self.r = self.full.shape[0]
        self.active = [t for t in st.S if t != 0 and v2i(abs(t)) < u]
        self.tgt = [x for x in self.full[:, st.S.index(0) * (N + 1)]]
    def data(self, offs):
        """data matrix with level-i columns multiplied by 2^{offs[i]}."""
        st = self.st; N = self.N; blocks = []
        for t in self.active:
            i = v2i(abs(t)); s = st.S.index(t) * (N + 1)
            blk = self.full[:, s:s + N + 1]
            if offs[i]: blk = np.array([[(x << offs[i]) & MASK for x in row] for row in blk], dtype=object)
            blocks.append(blk)
        return np.concatenate(blocks, axis=1)
    def need_additive(self, offs=None, p=1):
        """min P such that p_i = P - offs[i] is feasible."""
        offs = offs or [0] * self.u
        piv, left = smith_pivots(self.data(offs), self.tgt)
        return need_from(piv, left, p)
    def feasible(self, pvec):
        P = max(pvec); nd = self.need_additive([P - x for x in pvec])
        return nd is not None and nd <= P
    def need_balanced(self, hi=None):
        """min m such that p_i = m * 2^{u-1-i} (i.e. p_i proportional to D_i) is feasible."""
        u = self.u
        if u == 1: return self.need_additive()
        lo = 0; hi = hi or self.need_additive()        # uniform need is feasible => m=hi feasible
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if self.feasible([mid * 2 ** (u - 1 - i) for i in range(u)]): hi = mid
            else: lo = mid
        return hi

# ---------------- closed formulas
def formula_D9(a, u, M):
    """paper's taper h=ceil(3M/u): exact loss_i=-v2(E_t(0)); exact need of formula (D.9) (reviewer's count);
    paper's stated increment 8hD_t+8M(a+1)."""
    L = 2 ** a; h = -(-3 * M // u); q = [h * max(0, u - i) for i in range(a + 2)]
    out = {}
    for i in range(u):
        D = L // 2 ** i
        loss = D * sum(Fr(j, 2 ** j) * (q[i] - q[i + j]) for j in range(1, a - i + 1))
        out[i] = (int(loss), int(1 + loss + (q[i] - 1) * (a + 1 - i)), 1 + 8 * h * D + 8 * M * (a + 1))
    return out
def witness_g(a, u, M):
    """depth g_t of the product witness F_m (rev3), per level."""
    L = 2 ** a; N = M * (2 * L - 1)
    wsum = sum(Fr(L, 2 ** i) * (u - i) for i in range(u)); h = int(N // wsum)
    prof = [h * max(0, u - k) for k in range(a + 2)]
    return {i: int(Fr(L, 2 ** i) * sum(Fr(j, 2 ** j) * (prof[i] - prof[i + j]) for j in range(1, a - i + 1)) + prof[i] * (a + 1 - i)) for i in range(u)}
