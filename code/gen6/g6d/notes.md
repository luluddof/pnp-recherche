# G6d: per-point errors E_i in the OpenAI/O'Donnell interpolation scheme

## 1. Requirement (generalised Prop. 3.2 setting)
Single cube {0,1}^l, J(t) a 2^l x 2^l matrix polynomial, deg_t J <= d, J(0) = X^{(x)l}.
Points tau_1..tau_N (distinct, nonzero), point k has a direction c(k) in [l] (m_c points per direction).
J(tau_k) = K_k + E_k, K_k c(k)-nice. We only use the WEAK form: K_k[x,y] = 0 if x_{c(k)} != y_{c(k)}
(this is all Lemma 4.3 needs). A(t) = L (I_C (x) J(t)) R'; output error rank r_k <= |C| rank E_k = M rho_k,
rho_k := rank(E_k)/2^l.
Cost (paper accounting): sum_k 2^{O(n/q)} n^{O(M rho_k)} (n/h_k)^{O(h_k)},  with sum_k h_k > deg P (<= n d).
Exponent ~ max_k [ rho_k M log n + h_k log(n/h_k) ] + n/q, M >= n, l <= log n / (2 log q).

## 2. Lemma (subcube / determinant lemma) -- PROVED, general (non-commutative) J
For every set S of directions, |S| = s, with N_S = #points whose direction lies in S:
      sum_{k : c(k) in S} rank(E_k)  >=  (N_S - d) * 2^{l-s}.
Proof. A' = {x : x_S = 0}, B' = {y : y_S = 1}. Every (x,y) in A' x B' differs on all of S, so for c(k) in S,
K_k[A',B'] = 0 and E_k[A',B'] = J(tau_k)[A',B'] =: M(tau_k). M(0) = X^{(x)l}[A',B'] is a permutation matrix
(x -> complement), so det M(t) is a nonzero polynomial of degree <= d 2^{l-s}. nullity M(tau) <= ord_tau det M
(make m columns of M(tau) zero by a constant column operation; those columns of M(t) are divisible by t - tau).
Hence sum_k (2^{l-s} - rank M(tau_k)) <= d 2^{l-s}, and rank E_k >= rank M(tau_k).  QED
Corollaries.
 (a) one point per direction, S of size d+1:  max_k rank E_k >= 2^{l-d-1}/(d+1), i.e. rho >= 2^{-d-1}/(d+1).
 (b) S = {c}: m_c <= d/(1 - 2 rho). Averaging over all s-sets with 2^s = 1/(2 rho):
        N <= 2 d l / log2(1/(2 rho)).
     Hence h >= n d / N >= n log2(1/(2 rho)) / (2 l)   -- d cancels: more degree buys only
     proportionally more points.
 (c) restriction monotonicity: averaged normalised error rho_bar(l,d) >= rho_bar(s,d) for s <= l
     (disjoint subcube blocks A'_a x B'_a, ranks add).
 (d) G4c's common-E Hamming-ball bound is NOT implied and does NOT hold for per-point errors (see 3).

## 3. Constructions / numerics (commutative = Hadamard-diagonal class)
Reduction (exact): J(t) = H diag(g_z(t)) H; with g_z = chi(z) - t chi(z) r_z(t), r_z in P_{d-1}:
   rank E_k = #edges (z,z^k) with (r_z + r_{z^k})(tau_k) != 2/tau_k      (chi disappears).
 * Explicit degree-1 family, all l (d1check.py, exact incl. 2^l x 2^l matrix ranks for l<=6):
   r_w = 1/tau_1 + (-1)^{w_1} (-1)^{|w_{-1}|} sum_{j>=2} w_j (2/tau_j - 2/tau_1).
   Direction 1 exact, every other rank E_k = 2^{l-2} (rho = 1/4) for l = 2..10.
   Common-E at d = 1 needs rho >= ball(l,(l-2)/2)/2^l -> 1/2. So per-point errors strictly help.
 * MILP (HiGHS; optimal only under |values|<=30; bad set re-verified over Q), total bad edges T(l,d):
   solved to optimality: (2,1)=1 (3,1)=4 (3,2)=1 (4,1)=12 (4,2)=4 (4,3)=1 (5,3)=4 (5,4)=1,
   all = (l-d)2^{l-d-1} (d exact directions + rho = 2^{-d-1} on the rest).
   BUT (5,2): best found 14 (dual bound 8, 1200 s); with directions 1,2 forced exact the dual bound is 13 > 12,
   so the pattern "d exact directions + one bad copy per edge" BREAKS at l-d = 3, d = 2 (also predicted by an
   equation count). (6,2): best found 44 (frac 0.23). Large-l behaviour for d >= 2 is NOT determined.
   rigorous lower bounds come only from the Lemma (e.g. T(5,2) >= 7, T(6,2) >= 16).

## 4. Running time
Lower bound for the whole framework (paper's cost accounting, deg P bound n d used, M >= n, l = O(log n/log log n)):
   need rho M log n <~ target => log(1/rho) >= (2 - o(1)) log log n, so h >= n log log n / l
   => exponent >= Omega(n (log log n)^3 / log n).   (independent of d and of #points per direction)
IF some construction had rho = 2^{-Theta(d)} uniformly in l (one point per direction): d ~ C log log n,
   time 2^{O(n (log log n)^3 / log n)} -- matching the lower bound, vs OpenAI's 2^{O(n (log log n)^2/sqrt(log n))}.
   More generally, d(l, rho) = o(sqrt(l log l)) at rho = l^{-C} would beat sqrt(log n). OPEN.
Unconditionally (from what is proved): the sqrt(log n) loss is NOT forced by the per-point framework
(the G4c lemma does not extend); the framework cannot go below 2^{n/log n}.

## 5. Self-review / unproved points
- No construction with rho = 2^{-Theta(d)} for l >> d is known; the tempting formula (l-d)2^{l-d-1} fails
  (structurally) at (5,2). Linear "parity-split" ansatzes provably stall at
  rho = 1/4 (they need one instance vanishing at many extra points, impossible: each g_z has <= d roots).
  It is quite possible that rho(l,d) -> const or -> 1/2 for fixed d >= 2 as l -> infinity.
- Lower bound uses only weak niceness and the given L(I (x) J)R' shape; it treats output rank as |C| rank E_k
  (true cost could be lower if L, R' collapse rank) and assumes the interpolation uses deg P <= n d.
- MILP optimality relies on variable bounds |values| <= 30 (not a proof of a lower bound); lower bounds are
  only the lemma. Upper bounds are exact (rational re-verification).
