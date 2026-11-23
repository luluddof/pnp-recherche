"""Candidate complexity measures computed from a truth table (list/array of 0/1, length 2^n, input j = bit j of row)."""
import numpy as np, itertools, math

def tt_arr(f, n):
    return ((f >> np.arange(1 << n, dtype=object)) & 1).astype(np.int8) if n > 6 else \
        ((np.int64(f) >> np.arange(1 << n)) & 1).astype(np.int8)

def walsh(tt):
    a = (1 - 2 * np.asarray(tt, dtype=np.float64)).copy()
    h = 1
    while h < len(a):
        a = a.reshape(-1, 2, h)
        a = np.stack([a[:, 0] + a[:, 1], a[:, 0] - a[:, 1]], 1).reshape(-1)
        h *= 2
    return a / len(a)

def anf_degree(tt):
    a = np.asarray(tt, dtype=np.int8).copy(); h = 1
    while h < len(a):
        a = a.reshape(-1, 2, h); a[:, 1] ^= a[:, 0]; a = a.reshape(-1); h *= 2
    idx = np.nonzero(a)[0]
    return max((bin(i).count('1') for i in idx), default=0)

def spectral_L1(tt):  # log2 of Fourier L1 norm
    return math.log2(np.abs(walsh(tt)).sum())

def fourier_entropy(tt):
    w = walsh(tt) ** 2; w = w[w > 1e-15]
    return float(-(w * np.log2(w)).sum())

def total_influence(tt):
    t = np.asarray(tt); n = int(math.log2(len(t))); idx = np.arange(len(t))
    return sum(float((t != t[idx ^ (1 << j)]).mean()) for j in range(n))

def max_sensitivity(tt):
    t = np.asarray(tt); n = int(math.log2(len(t))); idx = np.arange(len(t))
    return int(sum((t != t[idx ^ (1 << j)]).astype(int) for j in range(n)).max())

def gf2_rank(M):
    M = M.copy() % 2; r = 0; rows, cols = M.shape
    for c in range(cols):
        piv = np.nonzero(M[r:, c])[0]
        if len(piv) == 0: continue
        p = r + piv[0]; M[[r, p]] = M[[p, r]]
        nz = np.nonzero(M[:, c])[0]; nz = nz[nz != r]
        M[nz] ^= M[r]; r += 1
        if r == rows: break
    return r

def comm_matrix(tt, A):
    t = np.asarray(tt); n = int(math.log2(len(t))); B = [j for j in range(n) if j not in A]
    M = np.zeros((1 << len(A), 1 << len(B)), dtype=np.int8)
    for x in range(len(t)):
        a = sum(((x >> j) & 1) << i for i, j in enumerate(A)); b = sum(((x >> j) & 1) << i for i, j in enumerate(B))
        M[a, b] = t[x]
    return M

def max_partition_rank(tt):  # max over balanced partitions of log2 GF(2)-rank of the communication matrix
    n = int(math.log2(len(tt)))
    return max(math.log2(max(1, gf2_rank(comm_matrix(tt, list(A))))) for A in itertools.combinations(range(n), n // 2))

def subfunction_counts(tt):
    """For each free set S (|S|>=1), number of distinct subfunctions on S under the 2^{n-|S|} assignments to the rest."""
    t = np.asarray(tt); n = int(math.log2(len(t))); out = {}
    for k in range(1, n + 1):
        for S in itertools.combinations(range(n), k):
            R = [j for j in range(n) if j not in S]
            subs = set()
            for r in range(1 << len(R)):
                base = sum(((r >> i) & 1) << j for i, j in enumerate(R))
                subs.add(tuple(int(t[base + sum(((s >> i) & 1) << j for i, j in enumerate(S))]) for s in range(1 << k)))
            out[S] = len(subs)
    return out

def nechiporuk_like(tt):
    """Sound-for-circuits variant: max over S of log2(#subfunctions on S). (A size-s circuit has at most
    (c s)^{2s} restrictions, so s log s >= Omega(this); this quantity is <= n - |S| <= n, hence useless beyond n.)"""
    c = subfunction_counts(tt)
    return max(math.log2(v) for v in c.values())

def nechiporuk_sum(tt):
    """Formula-style Nechiporuk sum over a best partition into singletons+pairs -- NOT sound for circuits."""
    c = subfunction_counts(tt); n = int(math.log2(len(tt)))
    return sum(math.log2(c[(j,)]) for j in range(n))

def onset_components(tt):  # topology: # connected components of on-set + off-set in the hypercube graph
    t = np.asarray(tt); n = int(math.log2(len(t))); seen = np.zeros(len(t), bool); comps = 0
    for v in range(len(t)):
        if seen[v]: continue
        comps += 1; stack = [v]; seen[v] = True
        while stack:
            u = stack.pop()
            for j in range(n):
                w = u ^ (1 << j)
                if not seen[w] and t[w] == t[u]:
                    seen[w] = True; stack.append(w)
    return comps

def essential(tt):
    t = np.asarray(tt); n = int(math.log2(len(t))); idx = np.arange(len(t))
    return sum(int((t != t[idx ^ (1 << j)]).any()) for j in range(n))

def xor_cover(tt):
    """NEW candidate ('XOR-defect'): F2-degree-weighted defect from affine functions.
    D(f) = (#essential vars - 1) + sum_{d>=2} (#distinct degree-d ANF monomials' supports needed, compressed)
    Concretely: number of nonlinear ANF monomials after taking the minimum over all
    invertible linear changes of the *output-side* is too costly; we use a cheap proxy:
    D(f) = (ess-1) + log2(1 + #nonlinear ANF monomials) * (deg-1)."""
    a = np.asarray(tt, dtype=np.int8).copy(); h = 1
    while h < len(a):
        a = a.reshape(-1, 2, h); a[:, 1] ^= a[:, 0]; a = a.reshape(-1); h *= 2
    idx = np.nonzero(a)[0]
    nl = sum(1 for i in idx if bin(i).count('1') >= 2)
    deg = max((bin(i).count('1') for i in idx), default=0)
    return (essential(tt) - 1) + math.log2(1 + nl) * max(0, deg - 1)

ALL = dict(ess_minus1=lambda t: essential(t) - 1, anf_deg=anf_degree, log_L1=spectral_L1, fourier_H=fourier_entropy,
           tot_inf=total_influence, max_sens=max_sensitivity, part_logrank=max_partition_rank,
           sub_max=nechiporuk_like, nech_sum=nechiporuk_sum, components=onset_components, xor_defect=xor_cover)
