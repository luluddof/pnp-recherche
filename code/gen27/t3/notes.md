# G27-T3 notes (start 08:08 UTC). Corner problem F(N,q), delta(eta).
Tags [P] proof sketch (unreviewed), [E] exact computation, [N] numeric, [H] heuristic.

## 0. Immediate remarks (08:10-08:35)
[P] PLAIN Pauli truncation is useless at the corner: K = (1/N) sum_{P in S} P (x) P^dag has ||K|| <= q/N, SWAP unitary,
    so for q < N SWAP-K is invertible: p = N^2 exactly (any N, any S, qubit or qudit Paulis). Coefficients MUST be
    re-optimised (consistent with B2 norm lemma Q >= N^delta).
[P] Qubit Pauli / Bell model (N = 2^k): P (x) P (P = X^a Z^b real) are simultaneously diagonal in the Bell basis
    Phi_{c,d}, eigenvalue (-1)^{a.d+b.c}; SWAP eigenvalue (-1)^{c.d}. For K = sum_{(a,b) in S} kappa_ab P(x)P:
    rank(SWAP-K) = #{v : phi(v) != beta(v)}, beta(c,d) = (-1)^{c.d} (bent), phi = sum_S kappa chi_S (q-Fourier-sparse).
    Matrix form: Phi = H C H (H Sylvester Hadamard N x N, C q-sparse), p = #{Phi != H}. Self-dual: (p,q) <-> (q,p).
[E] bell_check.py: Bell formula == explicit rank(SWAP-K) on random S, integer kappa, N=4 and 8 (all agree).
[E] bell_n4.py: EXACT Bell-model frontier at N=4 (all 69 ASp(4,2)-orbits of S, all info sets):
    q:1-3 -> p=6 ; q=4,5 -> 4 ; q=6..15 -> 1 ; q=16 -> 0.  Self-dual points (1,6),(4,4),(6,1).
    (General model is better at q=3: known (p,q)=(5,3)@4, so Bell is a strict submodel already at N=4.)

## 1. Bell "polynomial method" family (bell_poly.py) [E]
phi(x,y) = f(m), m = x.y in {0..k} (integer), f = Lagrange interpolant of (-1)^m on a window A of d+1 integers.
q = sum_{u: F(u)!=0} C(k,u) 3^u, F(u) = sum_j c_j 4^{-j} C(k-u,j-u) (c_j Newton coefficients of f);
p = sum_{m: f(m)!=(-1)^m} C(k,m) 3^{k-m}.  Formula == direct WHT count for all k<=7, all d (exact).
k=7 (N=128): (q,p) = (1,8256),(22,6178),(190,3343),(1156,1156),(3991,211),(9094,22),(8128,1)
  -> balanced p=q=1156 = N^1.454 at N=128 (cf. N^2/8 = 2048 Hurwitz-blocks; App B trivial there).
[E] bell_verify.py: explicit integer matrices n*D*(SWAP-K), rank mod 2147483629 and 2147483587 == Bell count for
    EVERY d at k=3,4,5 (N=8,16,32; e.g. N=32: (q,p) = (1,528),(16,349),(106,106),(376,16),(496,1)). verify_k5.txt.
    Schmidt rank of K = q exactly (P(x)P linearly independent).
[E] bell_bal.py: balanced max(p,q) of the family: k odd gives EXACT self-dual points p=q (d=(k-1)/2):
    N=32: 106 (N^2/8=128) | N=128: 1156 (2048) | N=512: 12826 (32768) | N=2^11: 143980 | N=2^13: 1630084 | N=2^15: 18575560
    exponents 1.35,1.45,1.52,1.56,1.59,1.61 ; N=2^24: 1.65 ; asymptotically -> 1/2+... = H(1/2)+log2(3)/2 = 1.792 (bad).
    => side result: BELOW the direct sum N^2/8 for every N = 2^k, 5 <= k <= 24 (JOURNAL: "no ansatz passes N^2/8 at any
       size" was about G23a's group ansatz). Finite-size only; asymptotically worse than App B.
[N] bell_asym.py: corner asymptotics of the family: eta = H(tau)+tau log2 3, delta = best two-sided KL window tail (bits):
    tau=.01: eta=.097 delta=1e-4 ; tau=.05: eta=.366 delta=.0024 ; tau=.1: eta=.63 delta=.0096 ; tau=.2: eta=1.04 delta=.039
    delta/eta^2 -> 0 like 1/log^2(1/eta): delta ~ eta^2/(C log^2(1/eta)) -- WORSE than App B (eta^2/4) at the corner.
    Reason: Gaussian (binomial) concentration of m = x.y + interpolating parity on a window of width w costs degree w.
[E] plain_trunc.py: plain truncation (coefficients 1/N) exhaustive N=4, greedy N=8: p = N^2 for every q < N, and p+q >~ N^2
    everywhere (N=8: q=32 -> p=44, q=56 -> p=20). Useless, as the norm argument says.
[N] bell_norm.py: Pauli 1-norm Q = sum|kappa| of the family: k=7: d=1,2,3 -> Q = 11, 49, 209 vs N^delta = 2.7, 4.9, 14.2;
    B2's norm lemma p >= (N^2-NQ)/(1+Q) holds (and is far from tight). Norms grow ~ q, as the lemma forces.

## 2. Structure results in the Bell model [P = proof sketch, unreviewed]
(B1) Bell <= Hadamard rigidity: Phi = H C H^T with C q-sparse has rank <= q, so p_Bell(q) >= R_H(q) (Sylvester H_N,
     rigidity = min #changed entries to reach rank <= q). Known: R_H(r) >= N^2/(4r) (Kashin-Razborov/de Wolf) -> only the line.
     Conversely a Bell linear corner would give R_H(N^eta) <= N^{2-c eta}: low-rank NON-rigidity of Hadamard with linear
     exponent (not known; the polynomial method -- Alman-Williams style -- gives only N^{2 - eta^2/log^2}).
(B2) Restriction monotonicity (Bell only): restrict phi to a random coset of a random symplectically-nondegenerate subspace
     of dim 2k' < 2k: beta stays bent (quadratic, Arf type irrelevant up to a character twist), sparsity <= q, average error
     FRACTION unchanged. So (eta, delta) achievable at N  =>  (lambda eta, lambda delta) at N^{1/lambda}, lambda = k/k' >= 1.
     => delta_Bell(eta)/eta is nondecreasing; the Bell corner slope = inf_eta delta/eta. (Fails in the general model:
     compressing R to W(x)W does not divide its rank by (N/N')^2 -- rank is not an error count. This is exactly why the
     general corner can be better than any finite-size ratio suggests, and why Bell data cannot be lifted.)
(B3) Row uncertainty: row x of phi - beta has Fourier support in pi_2(S) u {x}; if nonzero it has >= N/(q+1) errors, rows
     with x notin pi_2(S) are never exact => p >= (N-q) N/(q+1). Same for every Lagrangian splitting. Gives the line only;
     "row-tight" designs (errors = a coset per row) force q >~ N^2/u (computed), so tightness is impossible jointly.
(B4) THEOREM S (S_k-symmetric Bell constructions have a quadratic corner). N=2^k, phi invariant under simultaneous
     permutation of coordinates of x and y, q < min(C(k,u+1), C(k,k/4)). Then p >= c N^2 k^{-3/2} exp(-C u^2/k).
     Hence delta_sym(eta) <= C (eta/log(1/eta))^2 + O(log k/k): the polynomial family is optimal up to constants
     among symmetric Bell constructions; no linear corner.
     Proof: an S_k-orbit of characters whose type composition (u00,u01,u10,u11) has k - max_t u_t = dist >= u+1 has size
     >= min(C(k,u+1),C(k,k/4)) > q, so every character of S is within u of a pure type t^k (t in {00,01,10,11});
     phi*beta = sum_t sigma_t P_t with P_t polynomials of degree <= u in the type counts (n01,n10,n11) and
     sigma_t in {(-1)^{n11}, (-1)^{n10}, (-1)^{n01}, (-1)^{k+n00}} (= the 4 symmetric characters times beta).
     The agreement set Z is a union of type classes. If Z contained a cube of side 2u+2 in type space, each of the 8 parity
     classes c of (n01,n10,n11) would contain a (u+1)^3 grid, so sum_t eps_t(c) P_t == 1 identically for all 8 classes
     (tensor-grid unisolvence); but c -> eps_t(c) are 4 NON-trivial characters of (Z_2)^3, so summing over c gives 0 = 8.
     Hence Z misses a type within l_inf distance 2u+2 of the centre (k/4,k/4,k/4); that class alone has mass
     >= c k^{-3/2} exp(-C u^2/k) N^2.  QED (sketch).
     Moral: symmetric statistic + CLT => errors exp(-w^2/k) for cost C(k,w): quadratic corner (same shape as App B, whose
     cost is p(w) ~ exp(c sqrt w)). A linear corner needs an error mechanism with EXPONENTIAL (not Gaussian) tails at cost
     2^{O(w)}: not available in any symmetric/graded family I know.

## 3. Clean operator form of the family (08:40) [E/P]
Per qubit pair SWAP_2 = 1 - 2E_i (E_i singlet projector); the E_i commute, M := sum_i E_i (singlet number),
SWAP_N = prod (1-2E_i) = (-1)^M.  K = f(M), f of degree d interpolating (-1)^m on a window W of d+1 integers:
 rank(SWAP-K) = sum_{m notin W} C(k,m) 3^{k-m}  (Bell states with m singlets),  Schmidt rank(K) = #Pauli strings of
 weight <= d = sum_{u<=d} C(k,u) 3^u (all coefficients nonzero; checked by formula F(u) and direct WHT).
 Window [0,d] = truncation of prod(1-2E_i) at degree d. k odd, d=(k-1)/2: p = q = sum_{u<=d} C(k,u)3^u (self-dual).
 Qudit-site version (N = n^k): SWAP_n = 1 - 2A_i (A_i antisym. projector), same construction, q <= sum C(k,u)(n^2-1)^u,
 p exact; qudit_asym.py: delta/eta^2 ~ 0.01-0.02 at eta in [0.05,0.6] for n = 2,3,5,8 (all CLT-quadratic, App B ~0.25).
 bell_bal2.py: exponent of the self-dual point: k=25: 1.672, 31: 1.691, 41: 1.711, 61: 1.734 (-> 1.792).
Tensor/blocks (requested idea SWAP_N = SWAP_a (x) SWAP_b) [P]: (R_a+C_a)(x)(R_b+C_b) with R = R_a(x)SWAP_b + C_a(x)R_b,
 C = C_a(x)C_b gives delta = min(alpha delta_a, (1-alpha) delta_b), eta = alpha eta_a + (1-alpha) eta_b
 (N = N^alpha * N^{1-alpha}): a corner slope c becomes <= c/2; "C exact on one factor" (C_a = SWAP_a) gives
 T >= 2 alpha: moves away from the corner. Tensoring cannot create a linear corner (B2 agrees).
MILP (bell_milp.py, HiGHS, |kappa|<=4) at N=8 q=2..4: only finds p=28; LP bounds useless (gap 89%). Abandoned.

## 4. Lower-bound side (B) -- what I could and could not do
[P] Line, short proof (any field): V = ker(1-Theta), Theta(M) = sum_l A_l M^T B_l^T. For M in V with col(M) in U (dim a),
    row(M) = row(Theta M) lies in sum_l B_l U (dim <= qa). So aN - p <= dim(V n U(x)F^N) <= q a^2 for every a;
    a = N/(2q) gives p >= N^2/(4q). (Same constant class as Theorem C; it is a rank/dimension argument, so by Theorem B it
    cannot be pushed past the line, and a linear corner p = N^{2-c eta}, q = N^eta lies ABOVE the line (pq = N^{2+(1-c)eta}).)
=> Excluding a linear corner (delta <= C eta^2, or just delta = o(eta)) needs a beyond-line bound at q = N^eta: impossible
   for rank methods (Theorem B); spectral/norm methods (Theta_s lemma + odd traces) are capped at N^2 log N / Q and only
   force Q >= N^delta (B2 lemma; the Bell family has Q ~ q >> N^delta, consistent) -- they never see q.
[P] Bell: p_Bell >= R_H(q) (Hadamard rigidity); best known R_H(r) >= N^2/(4r): line only. Proving delta_Bell <= C eta^2
    would need a new low-rank rigidity-type bound for Fourier-sparse approximants of a bent function (open as far as I know).
[H] Norm-tail tradeoff (Bell): |E[beta phi]| <= ||phi^||_1 / N, so if ||phi^||_1 = A << N then E[|phi| 1_bad] >= 1-eps-A/N:
    phi must be ~ N^delta on the bad set. Consequence: "commuting-reflection" constructions K = f(sum_j g_j) with g_j
    low-complexity Boolean, deg f = w, A <= C(K,w) 2^{O(w)}, CANNOT have super-exponential (Poisson) tails of m = sum g_j at
    small eta; tails must decay no faster than |f| grows outside the window. This kills the most naive "exponential-tail"
    route to a near-linear corner inside Bell (the would-be slope was ~ c/log(K/w) anyway).
NOT DONE: qudit Weyl subsets for odd prime powers (P(x)P^dag do not commute, no Bell reduction; plain truncation is useless
 by the norm argument for every N; the commuting qudit-site family of sect. 3 is the natural replacement and is quadratic).

## VERDICT (08:41, stopped early on coordinator deadline)
No linear corner found. Every family I built or analysed is quadratic up to logs (Bell/singlet-number family
delta ~ eta^2/(C log^2(1/eta)), worse than App B's eta^2/4; proved for S_k-symmetric Bell). No proof of delta <= C eta^2:
it needs a beyond-line bound, which rank methods (Thm B) and norm/spectral methods (B2/G26) cannot give.
Corner route: OPEN, and as hard as the line. Side result: finite-size self-dual points below N^2/8 for N = 2^5..2^24.
Unfinished: Theorem S for general (non-Bell) S_k-symmetric operators; exact Bell optimum at N=8,16; qudit Weyl subsets.
