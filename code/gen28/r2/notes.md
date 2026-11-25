# G28-R2 review of G28-C2 (start 12:17 UTC)
## 12:25 hand check of lemmas
- L1: K inv under pi(x)pi <=> C(pi a,pi b)=C(a,b) <=> C intertwiner; im C, im C^T S_k-stable; any rank decomposition has
  A_l in im C (so EXISTENCE of a decomposition with bounded-degree factors; enough since rank(SWAP-K) depends on K only). OK.
- L2: = Dvir/Clausen-Meier type fact (g_{mu nu nu'}>0 => |nu cap nu'| >= mu_1). Proof via M^mu < Ind(triv x C S_j) + LR OK.
- L3: correct, BUT NOT SHARP: threshold 2u can be replaced by u. Y={f(t)-f(t')>u}: K(Y) in {f(s)-f(s') > -u}, X={f(s)-f(s')<-u}.
  => p >= N^2 Pr[L1-L2>u] ; banded: (N-w)(N-w-1)/2 (direct: fixed permutation submatrix in min-rank completion).
  => asymptotic constant improves by 4: delta_sym <= (d^2 ln d/(2(d-1))) eta^2/ln^2(1/eta).
## 12:45 numerics
- exact_points.py (own build: K = sum_{nu!=nu'} P_nu(x)P_nu' + sum_nu SWAP.M_nu, M_nu = (f/k!)^2 sum chi_nu(sig pi) pi(x)sig,
  integer chars of two-row shapes, mod 1000003 and 998244353): (20,20)@8, (101,35)@16, (452,56)@32 exact; float SVD same
  (gap 2 vs 1e-13), K is pi(x)pi-invariant. Dominated: (62,12)@16 [G27-B1], (32,32)@16 HR tiling [R1], (349,16)@32 Bell [T3/R1].
- lemma_check.py (d,k)=(2,4),(3,3),(2,5),(4,3 partial): conjugation-isotypic mu random ops: Lemma 2 (common alpha |-k-j, band)
  0 violations. Lemma 3: leak P_X K P_Y = 0 at threshold u (improved), nonzero at u-1; C2's 2u gives empty Y at these k.
- asym.py: finite k reproduces C2 (k=1600, eta=.1: .0039). k->oo limit / formula ratio 0.19 (eta=.1) .. 0.58 (1e-6) -> formula
  is a valid (1+o(1)) upper bound, slow convergence. Improved threshold: constant 5.545 -> 1.386 (d=2).
- Bell span-r corollary: Witt puts all factors on the SAME s+t<=r qubits -> K = K_r (x) 1, so p >= M_- 4^r = (N^2 - N 2^r)/2
  (eigenspaces of SWAP on remaining qubits). Much stronger than C2's N^2 Pr[L1-L2>2r].
- band_tight.py: banded min rank = (N-w)(N-w-1)/2 EXACTLY (explicit K: cancel |i-k|<=w swaps, rank-1 2x2 blocks on the rest),
  N=5..10 all w. C2's (N-2w)(N-2w-1)/2 is a valid but loose bound, and is wrong as written for 2w > N (needs (.)_+).
- bell_span.py: optimum p over random span-r Pauli sets >= (N^2-N2^r)/2 (k=3,4,5); C2's bound is 0 at these sizes.
- weyl_mine.py: d=3 k=1,2 (two different anti-symplectic involutions), d=5 k=1: SWAP monomial in common eigenbasis, N fixed
  points, p_opt = (N^2-N)/2. Claim 4 confirmed.
- hookmin_mine.py: two-row minimiser (k<=28, j<=k/2) and hook lemma f^mu sqrt(e) >= C(k,j) (j<=k/4): 0 violations.
  eta_0(d) branching argument checked by hand (d=2: .270, d=3: .076).
- R2 by hand: phi(c,d) = sum_{(a,b) in S} phihat (-1)^{a.c+b.d} rank <= q, differs from H (N x N) in p entries => R_H(q) <= p.
  Direction correct. Random translation (beta(c+a,d+b) = beta (-1)^{a.d+b.c+a.b}) gives prob-rank statement: correct.
## 13:00 C3 (Theorem S'') proof idea
- osr = sum_{c_U!=0}(n^2-1)^{|U|} (adapted basis 1 + sl_n per site, blocks by U) and p = N^2 mu_b(g != parity): c3_check*.py
  n=2 k=3,4; n=3 k=2,3 exact match (rank, realignment rank).
- mixture of subcubes (fix 0 w.p. 1-2b, free w.p. 2b) has marginal b: OK. Smolensky over R on agreement set: OK.
- constant: leading order only. Exact k->oo bound (optimised over m/k) = 0.574 eta^2 (1+O(eta^2)): ratio 1.006 at .3, 1.017 at .5
  (c3_asym.py). Model is a submodel of C2's u-local corollary (u = deg g) which gives n^2 ln n/(2(n-1)ln^2(n^2-1)) eta^2
  (improved threshold): C3's subcube argument is sharper by factor n.
