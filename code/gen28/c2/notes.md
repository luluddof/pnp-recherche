# G28-C2 notes (start 09:44 UTC). Corner problem, Theorem S beyond Bell.
Tags: [P] proof (unreviewed, written carefully), [S] sketch, [E] exact computation, [N] numeric, [H] heuristic.

## 09:44-10:05 reading + main idea
Key idea (found 10:00): do NOT diagonalise. Use the FIRST-FACTOR / SECOND-FACTOR Schur-Weyl sectors (GL_d x S_k on each
(C^d)^{(x)k}) and a 1-Lipschitz potential on sector labels (first row nu_1). SWAP exchanges sector labels (nu,nu')->(nu',nu);
low-"S_k-degree" factors move nu_1 by <= u. A compression P_X(SWAP-K)P_Y = P_X SWAP P_Y then has full rank dim Y.

## 10:05 Literature check (coordinator lead: low-rank Hadamard rigidity) [L = web search snippets, papers NOT read; arxiv.org blocked]
- Alman-Williams 2017 (arXiv 1611.05558): H non-rigid at rank N^{1-c eps^2} with N^{O(eps log 1/eps)} changes/row; also "non-trivial
  rigidity upper bounds for smaller target rank" (exact statement not retrieved).
- Alman-Liang 2025 (arXiv 2502.19580) "Low rank matrix rigidity: tight lower bounds and hardness amplification": over finite
  fields R_H(c1 log N) >= ~N^2/2 (tight); lower bounds at rank log^{1+eps} N would amplify to much stronger bounds; low-rank
  rigidity UPPER-bound improvements are tied to constant-depth linear circuits for H (Alman's depth-2 N^{1.476}).
- Alman 2026 (arXiv 2608.06592): over F_3, changing 1% of entries of H cannot reduce rank to <= log^2 N/80: "first constant-
  fraction rigidity lower bound for an explicit family at superlogarithmic rank over any field".
  => state of the art for EXCLUDING low-rank approximants of H with a constant error fraction is rank ~ log^2 N (finite field);
     over R the explicit bound is de Wolf/Kashin-Razborov N^2/(4r). Our corner needs error fraction N^{-delta} at rank N^eta.

## 10:20 THEOREM S' (all S_k-invariant K; proof) [P] + checks [E/N]
Setting. N = d^k, H = (C^d)^{(x)k} (first factor) (x) (C^d)^{(x)k} (second). S_k acts on sites; K is site-permutation invariant:
(pi (x) pi) K (pi (x) pi)^{-1} = K. osr(K) <= q. p = rank(SWAP - K).
Schur-Weyl on ONE factor: (C^d)^{(x)k} = (+)_nu S_nu (x) W_nu (nu |- k, <= d rows), Pi_nu = sector projector.
Pr_SW(nu) = f^nu dim W_nu / d^k (law of the RSK shape of a uniform word in [d]^k); L(nu) := nu_1.
Def. An operator A on (C^d)^{(x)k} has S_k-DEGREE <= u if it lies in the sum of the isotypic components, for the CONJUGATION
action A -> pi A pi^{-1}, of type mu with mu_1 >= k-u. Equivalently (L1) A in span{pi (X^{(x)(k-u)} (x) Y) pi^{-1}}; in
particular every sum of u-local operators has degree <= u.
LEMMA 1 (rank -> degree). If K is S_k-invariant with osr(K) <= q, then K = sum_{l<=q} A_l (x) B_l with all A_l, B_l of degree
 <= u(q) := max{k - mu_1 : mu |- k, l(mu) <= d^2, f^mu <= q}.
 Proof. Expand K = sum_{a,b} C(a,b) E_a (x) F_b in product operator bases (E_a = (x)_i E_{a_i}, a in [d^2]^k). S_k-invariance of K
 <=> C(pi a, pi b) = C(a,b) <=> C commutes with the permutation action on V^{(x)k}, V = C^{d^2} (labels). So im C and im C^T are
 S_k-stable, of dim = rank C = osr(K) <= q, hence contained in the sum of isotypics mu with f^mu <= q (an S_k-stable subspace
 meeting the mu-isotypic has dim >= f^mu). Take a rank decomposition C = sum alpha_l beta_l^T with alpha_l in im C, beta_l in
 im C^T; A_l = sum alpha_l(a) E_a. The label permutation action = conjugation by site permutations. QED
LEMMA 2 (degree -> band). If A has degree <= u then Pi_nu A Pi_nu' = 0 unless nu, nu' contain a common alpha |- k-j, j <= u;
 in particular |nu_1 - nu'_1| <= u.
 Proof. A -> Pi_nu A Pi_nu' is conjugation-equivariant, and Hom(S_nu'(x)W_nu', S_nu(x)W_nu) = (S_nu (x) S_nu'^*) (x) Hom(W,W)
 with S_k acting on the first factor. So the block is nonzero only if S_mu < S_nu (x) S_nu' (Kronecker) for some mu_1 = k-j >= k-u.
 S_mu < M^mu < Ind_{S_{k-j} x S_j}(triv (x) C[S_j]); tensor identity: S_mu (x) S_nu' < Ind_{S_{k-j}xS_j}((triv (x) C S_j) (x) Res S_nu')
 = sum over alpha subset nu' (|alpha| = k-j) of Ind(S_alpha (x) [C S_j]^{mult}); by Littlewood-Richardson every constituent
 S_nu contains alpha.  QED
LEMMA 3 (grading lemma, general). Let C^N = (+)_t H_t (orthogonal, same decomposition used on both tensor factors) and f a
 real function on labels. Suppose every A_l, B_l in K = sum A_l (x) B_l satisfies Pi_t A Pi_t' = 0 whenever |f(t)-f(t')| > u.
 Then rank(SWAP - K) >= sum_{f(t)-f(t') > 2u} dim H_t dim H_t'.   (NO bound on the number of terms is needed.)
 Proof. Blocks H_t (x) H_t'. Y := sum of blocks with f(t)-f(t') > 2u; X := SWAP(Y) = blocks with f(t')-f(t) > 2u.
 K maps block (t,t') into blocks (s,s') with f(s) >= f(t)-u, f(s') <= f(t')+u, so f(s)-f(s') > 0: never in X. Hence
 P_X (SWAP-K) P_Y = P_X SWAP P_Y, an isometry from Y onto X: rank >= dim Y.  QED
THEOREM S'. S_k-invariant K on (C^d)^{(x)k} (x) (C^d)^{(x)k}, osr(K) <= q:
     rank(SWAP - K) >= N^2 * Pr[ L_1 - L_2 > 2 u(q) ],   L_1, L_2 iid first-row lengths of SW-random shapes
 (= longest weakly increasing subsequence of a uniform random word in [d]^k).   [Lemmas 1+2+3 with f = nu_1.]
 Quantitative u(q) [P]: hook formula gives f^mu = C(k,j) f^{mubar} / E, E = prod_{c<=mu_2}(1+(mu'_c-1)/(mu_1-c+1)) <= e^{j/(k-2j+1)},
 so f^mu >= C(k,j)/sqrt(e) for j = k-mu_1 <= k/4; for mu_1 < 3k/4 (and l(mu) <= d^2, so mu_1 >= k/d^2) branching f^mu >= f^nu,
 nu = (mu_1, any |nu|-subdiagram of size floor(mu_1/3)), gives f^mu >= C(m, floor(mu_1/3))/sqrt e, m ~ 4k/(3d^2).
 Hence u(q) < u* := min{j: C(k,j+1) > sqrt(e) q} whenever this is < k/4 and q < 2^{c k/d^2}.
 Tail [P]: L >= max letter count, so Pr[L_1 >= k/d + t] >= Pr[Bin(k,1/d) >= k/d + t]; Greene + Doob: E L <= k/d + C_d sqrt k, so
 Pr[L_2 <= k/d + 2C_d sqrt k] >= 1/2. Hence p >= (1/2) N^2 Pr[Bin(k,1/d) >= k/d + 2u + 2C_d sqrt k].
 ASYMPTOTICS (fixed d, k -> oo, q = N^eta): u/k -> x, H_e(x) = eta ln d; delta(eta) <= 2 d^2 x^2/((d-1) ln d)
   ~ (2 d^2 ln d/(d-1)) eta^2 / ln^2(1/eta).   d=2: delta <= 5.55 eta^2/ln^2(1/eta) (asymptotically).
 => EVERY site-symmetric construction (not only Bell) has a quadratic corner (up to log^2). Proof is a pure sector/grading
    argument, no unisolvence, and it contains Theorem S (Bell) as a special case.
 Relation to Theorem G (G26-P3): u = 0 means every A_l, B_l commutes with S_k (term-wise symmetric): S' gives
    (N^2 - sum (f w)^2)/2, G gives (N^2 - sum f w^2)/2. So S' is a ROBUST Theorem G for Gamma = S_k ("approximately symmetric
    terms, degree <= u"), i.e. exactly the kind of inverse-type statement P3 asked for, in this group.
COROLLARY (locality, no symmetry needed) [P]: if all factors A_l, B_l are sums of u-local operators (any number q of terms),
 p >= N^2 Pr[L_1 - L_2 > 2u]. A construction with corner exponent delta needs factors of locality u >~ k sqrt(delta ln d)/C
 in EVERY tensor factorisation C^N = (C^d)^{(x)k}. Banded factors (grading H_t = C e_t, bandwidth w) give p >= (N-2w)(N-2w-1)/2.
CHECKS. band_check.py [N]: L1 dimension identity exact (k=4: 35,170,255,256 = predicted sums); L2 forbidden blocks <= 1e-12
 (k=4,5,6, all u). sprime_check.py [E]: bound holds on window family / truncations at N=4,8,16 (only u=0 non-vacuous there:
 N=16: bound 73, true optimum class >= 120). sprime_asym.py [N] (d=2): at k=1600, eta=.1: delta_max .0039 (bound) vs T3
 family delta ~ .001; eta=.3: .038 vs ~.01. Bound within factor ~3-5 of the best symmetric construction: the symmetric
 class has delta = Theta~(eta^2) (log factors aside). Finite-k: sqrt(k) offset dominates at small eta (k=100,400).

## 10:30 Task 1(b): diagonal ("Bell-type") models for qudits
(b0) [P] S_k-invariant diagonal models (any local SWAP-eigenbasis, any d) are a subclass of S' => quadratic corner.
(b1) [P] Weyl-Bell models for odd d: K in span{W_u (x) W_{Mu}}, M anti-symplectic involution (M(c,e)=(c,-e) <-> W (x) conj W).
  These commute (common eigenbasis = qudit Bell basis Phi_v = (W_v (x) 1)Phi), but SWAP Phi_{(c,e)} = phase * Phi_{(c,-e)}:
  SWAP is a monomial with N fixed points and (N^2-N)/2 two-cycles in that basis; each 2-cycle block [[-a,s],[s',-b]] has
  rank >= 1, so p >= (N^2-N)/2 for every K in the span (tight with q = N^2). weyl_check.py / weyl_opt.py [E]: d=3 k=1,2 and
  d=5 k=1: fixed points = N, p_opt = (N^2-N)/2 exactly. => Weyl-diagonal models are USELESS for odd d (Theorem-G level).
(b2) [S] A "true" Bell structure (SWAP in a commutative algebra spanned by d^2 linearly independent product operators X_a (x) Y_a)
  forces Y_a ~ X_a (each term commutes with SWAP) and X_a pairwise commuting or anticommuting; for invertible X_a this is an
  irreducible projective rep of an abelian group with +-1 commutator form => d = 2^m. So Bell (osr = Fourier sparsity,
  p = error count) exists only for qubit-like sites; for other d, diagonal models have osr not equal to a sparsity count.
(b3) For non-symmetric diagonal models nothing beyond the line is known; for d=2 it is the Bell model = Hadamard-rigidity-type
  problem (see task 2).

## 10:40 Task 2: general model at the corner -- what is known, reductions
KNOWN (upper bounds on delta = lower bounds on p):
 - delta <= eta (Theorem C, pq >= (2/9)N^2; dimension proof p >= N^2/(4q)) -- the line. [P+, Lean]
 - MC (signed monomial factors): p^2 q >~ N^4/8 => delta <= eta/2 there. Does NOT cover weighted monomials (Bell P(x)P with
   real coefficients) nor general K.
 - BN/norm: p >= c N^2 ln N/Q, forces Q >= N^delta; never sees q.
 - NEW (this note): S_k-invariant / u-local / banded classes: delta <= C (eta/log(1/eta))^2 (Theorem S', grading lemma).
NO argument giving delta <= (1-eps) eta in the general model is known to me, and none can be a class-B rank method (Thm B:
 delta <= (1-eps)eta means pq >= N^{2+eps eta}, beyond the line). The grading lemma IS a rank (compression) method; it beats
 the line only because the class restricts the factors. For general K one can build a bandwidth-1 grading from the joint
 Krylov spaces of {A_l, B_l, adjoints} (H_t = span of words of length t applied to a seed, orthogonalised), but dims grow like
 (4q)^t, mass sits in the top levels, and Lemma 3 gives only p >~ N^2/(4q)^3 (below the line). [P, easy]  The symmetric
 classes win because their sector statistic has CLT spread sqrt k with Gaussian tails while bandwidth u costs C(k,u) terms:
 that mismatch (Gaussian tail vs exponential cost) IS the quadratic corner. A linear corner needs a grading-free mechanism or
 a grading with EXPONENTIAL tails at exponential cost.

REDUCTIONS (N = 2^k; Bell model of G27-T3: phi q-Fourier-sparse on F_2^{2k}, p = #{phi != beta}, beta(c,d) = (-1)^{c.d} = H):
 (R1) [P] F(N,q) <= p_Bell(N,q): delta_gen >= delta_Bell. So "no linear corner in general" => "no linear corner in Bell".
 (R2) [P] Bell (p,q) => R_H(q) <= p (phi as N x N matrix in (c,d) is a sum of q rank-1 characters). Hence
      LINEAR BELL CORNER => R_H(N^eta) <= N^{2 - c eta} for all small eta: Hadamard NON-rigidity at polynomially small rank with
      a linear exponent. Stronger: by the translation symmetry beta(c+a,d+b) = beta(c,d)(-1)^{a.d+b.c+a.b} (shift+modulation keeps
      Fourier sparsity), a Bell pair (p,q) gives a PROBABILISTIC Fourier-q-sparse (hence rank-q) representation of H with per-entry
      error p/N^2: prob-rank_eps(H) <= q, eps = p/N^2.
 (R3) Best known low-rank NON-rigidity of H (upper bounds on R_H(N^eta)): polynomial method = T3's window family:
      delta ~ eta^2/(C log^2(1/eta)) [E, T3]. Alman-Williams 2017 mention "non-trivial upper bounds for smaller target rank"
      (statement not retrieved). Nothing linear is known to me. Moreover the two standard mechanisms are provably quadratic:
      (i) symmetric constructions (Theorem S'), (ii) real probabilistic polynomials in the AND variables z_i = c_i d_i:
      probabilistic degree of parity over R is ~ sqrt(k log(1/eps)) (Srinivasan-Tripathi-Venkitesh 2019/21 [L, recalled,
      near-tight for all symmetric functions]) => sparsity sum_{j<=w} C(k,j)4^j with w^2 ~ k log(1/eps): quadratic again.
      => A Bell corner (hence possibly step 2) would be a NEW non-rigidity theorem for Walsh-Hadamard in the low-rank regime,
         beating both known mechanisms. Not a known open problem in exactly this form, but adjacent to Alman-Liang 2025
         (low-rank rigidity of H is tied to constant-depth linear circuits for H; hardness amplification).
 (R4) Transfer of hypothetical non-rigidity: R_H(N^eta) <= N^{2-c eta} does NOT give a Bell corner (Bell needs the approximant to be
      Fourier-q-sparse, rank is not enough: a rank-1 phi = a(c)b(d) has Fourier sparsity |supp a^| |supp b^|), and does not give
      a general-model corner (Bell-diagonal K has osr = Fourier sparsity exactly, basis independent). So rigidity upper bounds do
      not transfer; only Fourier-SPARSE (sparse+sparse) decompositions of IP do. [P]
 (R5) Lower-bound side: "delta_gen = o(eta)" => "delta_Bell = o(eta)" <= "R_H(N^eta) >= N^{2-o(eta)}" (sufficient, not necessary).
      State of the art for excluding low-rank approximants of H: over R, R_H(r) >= N^2/(4r) (= line); over finite fields
      constant-fraction bounds only up to rank ~ log^2 N (Alman 2026, F_3; Alman-Liang 2025: rank c log N) [L]. Our need: error
      fraction N^{-o(eta)} at rank N^eta, i.e. near-constant-fraction rigidity at POLYNOMIAL rank: far beyond the frontier.
      So the rigidity route to "no corner" requires a breakthrough; the Fourier-sparse version (what Bell really needs) is a
      sparse-plus-Fourier-sparse decomposition question for IP for which only the uncertainty-principle line is known
      (row argument B3 of T3). I could not prove an equivalence with a named open problem; the clean statements are R2 and R5.

## 10:50 more on the coordinator's rigidity question [P/H]
- Is R_H(N^eta) <= N^{2-c eta} "easy"? Not that I can see. Natural constructions: block/Kronecker splits H = H_a (x) H_{k-a}
  (error 1/2), polynomial method (error e^{-w^2/k} at rank C(k,<=w)): both give delta = O(eta^2 / log^2). The latter is optimal
  among real probabilistic polynomials in the AND variables (prob. degree of parity) and among S_k-symmetric Bell
  constructions (Theorem S'); a linear low-rank non-rigidity would need a genuinely new idea (also: by Alman-Liang's
  connection, low-rank rigidity upper bounds for H interact with depth-2 linear circuits for H [L, not read]).
- Even if it were known, it would NOT transfer (R4): Bell needs Fourier-sparse approximants; the general model needs osr,
  which for Bell-diagonal K equals Fourier sparsity exactly.
- Converse (clean): Bell linear corner => (i) R_H(N^eta) <= N^{2-c eta}, (ii) prob-rank_{N^{-c eta}}(H) <= N^eta, both via
  explicit FOURIER-SPARSE approximants. So the Bell route to step 2 is literally "find a new non-rigidity theorem for
  Walsh-Hadamard at polynomially small rank with a linear exponent".
- Border: Theorem S' holds for border decompositions too (the band condition is a linear subspace condition on K, and the
  compression rank bound is uniform).

## 11:00 Task 3 numerics (small N)
[E] sector_constr.py: exact u=0 sector construction (Theorem-G optimum, factors commute with S_k; non-Bell):
    (p,q) = (20,20)@8, (101,35)@16, (452,56)@32 (p = (N^2 - sum f w^2)/2 exactly; q = dim of the S_k commutant = C(k+3,3)).
    Dominated by known points ((20,6)@8; Bell family (349,16)@32). Confirms S'/G at u=0 are tight (pred = exact).
[N] ap_sym.py (alternating projections Fl1-rank p / Fl2-rank q, optional S_k symmetrisation):
    N=4 general: recovers (6,1),(5,3),(4,4) (misses (1,6)): weak heuristic.
    N=16 S_4-symmetric, q=4: only p = 120 = N(N-1)/2 found (K ~ identity); p in 130..200 NOT converged (residual ~3e-2):
    AP + symmetrisation is unreliable (non-monotone). Not usable to estimate F at q ~ sqrt N.
[E] hook_check.py: f^mu sqrt(e)/C(k,j) >= 1.128 for all mu |- k <= 26 with j = k-mu_1 <= k/4 (hook lemma of S' confirmed).
[E] hook_min.py: for k <= 30 and j <= k/2, min{f^mu : mu_1 = k-j} is attained at the two-row shape (k-j,j) (all rows allowed).
    => u(q) = max{j : C(k,j)(k-2j+1)/(k-j+1) <= q} in practice. PROVED range of the asymptotic statement: eta < eta_0(d)
    ~ (4/(3d^2)) H_2(1/4)/log2 d (d=2: .27, d=3: .076) via the crude branching bound; beyond that the u(q) formula rests on
    the [N] two-row-minimiser fact.
[N] sprime_d3.py (d=3, N=3^k): k=400: eta=.05: delta<=.0041; .1: .0090; .2: .035; .3: .098 (delta/eta^2 ~ .9-1.6; T3 qutrit-site
    family has delta/eta^2 ~ .01-.02). Same quadratic shape.

## 11:10 Bell corollary of the locality theorem [P]
Clifford conjugation U (x) U preserves SWAP, osr and rank. If S (Pauli set of a Bell K) spans a subspace of F_2^{2k} of dim r,
a symplectic (Witt) change of basis maps it into the Paulis of r qubits, so all factors become r-local:
  p_Bell >= N^2 Pr[L_1 - L_2 > 2r].   With q <= 2^r this is quadratic (delta <= C r^2/k^2); so a Bell corner needs S spanning
  dim >~ k sqrt(delta) (no contradiction: q = N^eta >> k lets S span everything).
[N] sym_full_check.py: 400 random S_4-invariant K = sum_{r<=2} c_r (X_r (x) Y_r)^{(x)4} (osr <= 2): min p = 120 >= S' bound 73.
[N] ap_sym N=16 S_4-sym: q=8: p=180 found, 160..120 not (even the trivial 120): symmetrised AP is useless. General AP N=16
    q=4: p=112 FEASIBLE (residual 5e-10), p=119 not (non-monotone again).
[N] general AP N=16 q=4: p=112 feasible (res 5e-10); 104, 96 not (res > 3). => numeric point (112,4)@16:
    eta = .5, delta = 2 - log16(112) = .298: delta/eta = .60, delta/eta^2 = 1.19 -- one point cannot separate the laws.

## 11:20 Bell beyond-line attempt: column relaxation [P] + coset-tight counting [P]
PROP (column relaxation). For Bell (phi, S), T := pi_2(S) (|T| <= q):
   p_Bell >= sum_{c notin T} dist_Ham(chi_c, span_R{chi_t : t in T})  >=  Rcol(q) := min_{|T|<=q} sum_c dist(chi_c, span chi_T).
 (Row c of phi has spectrum in T; the row of beta is chi_c.)  Rcol is a "rigidity of H against approximants whose column space
 is spanned by q columns of H". Uncertainty gives dist >= N/(q+1) per column (the line).
LEMMA (coset-tight columns are rare). If column c is within m errors of span chi_T via an error function supported on a coset
 of a subspace W (|W| = m), then (c + W^perp) \ {c} subset T, s := |W^perp| = N/m, and c = t1+t2+t3 in >= (s-1)(s-2) ordered ways
 with t_i in T (t1 = c+u1, t2 = c+u2, t3 = c+u1+u2). Since sum_c r_3(c) <= q^3:  #{such c} <= q^3/((s-1)(s-2)) ~ q^3 m^2/N^2.
 For m = N/q^{1-eps}: at most ~q^{1+2 eps} << N columns. So the line can only be approached via NON-coset error patterns
 (|supp| x |spec| within a factor K = q^eps of N). An inverse uncertainty theorem with poly(K) losses (PFR-type, now known
 over F_2^n: Gowers-Green-Manners-Tao) might extend the count and give delta_Bell <= (1-eps) eta: beyond the line in Bell,
 BUT that would still not exclude a linear corner (step 2 needs only delta >= c eta for some c > 0; excluding it needs
 delta = o(eta), i.e. control of error patterns at slack K = q^{1-o(1)}, where no inverse theorem can help).
 Smolensky-type potentials (agreement set G => every chi_a|_G reduces to chi_{a+c+T}) give the Gaussian law for Hamming-ball T
 (agreement <= Pr[Bin(m,1/2) <= (m+w)/2]); for general T there is no potential: same obstruction as the grading lemma.
[E] rcol.py (exact circuit enumeration): N=16, q=3: Rcol = 100 ~ 0.39 N^2 for every T tried (d_c <= N/2 always via chi_c - chi_t;
    below N/2 only for c completing an affine subspace). Tiny q only; not informative about the corner (C(N,q) blow-up).

## 11:30 Heuristic principle from the grading lemma [H]
Corner exponent ~ (tail of the dimension profile of a grading) vs (cost of bandwidth):
 - product/CLT structure (sites, S_k sectors, Hamming weight): mass beyond distance u ~ e^{-u^2/k}, cost C(k,u) ~ e^{u log(k/u)}
   => delta ~ eta^2/log^2  (Theorem S', T3 families, Smolensky).
 - App B (Heisenberg Fock space, energy grading, dims p(n) ~ e^{c sqrt n}): sub-exponential growth -> quadratic (eta^2/4),
   consistent with T3's remark "cost p(w) ~ exp(c sqrt w)".
 - exponential-growth gradings (Krylov spaces of q generic operators: dims (4q)^t) do NOT obstruct: the lemma gives only
   p >~ N^2/(4q)^3. So if a linear corner exists, it should live on an exponential-growth ("free"/tree-like) structure --
   e.g. full Fock spaces over words, Cuntz/free-group-type operators -- not on tensor-product/site or bosonic structures.
   This is a constructive lead, not evidence.

## FINAL SUMMARY (G28-C2)
PROVED [P, unreviewed; lemmas numerically checked]:
 S'  (Theorem S beyond Bell). Every S_k-site-invariant K on (C^d)^{(x)k}(x)(C^d)^{(x)k} with osr <= q:
     rank(SWAP-K) >= N^2 Pr[L_1 - L_2 > 2u(q)], u(q) = max{k-mu_1 : f^mu <= q, l(mu) <= d^2} (~ two-row: C(k,u)(k-2u+1)/(k-u+1) <= q).
     => delta_sym(eta) <= (2d^2 ln d/(d-1)) eta^2/ln^2(1/eta) (1+o(1)), eta < eta_0(d). Holds for border decompositions.
     Ingredients: L1 rank->degree (im C is S_k-stable), L2 degree->band (Kronecker/LR), L3 grading lemma (compression).
 Grading lemma / locality corollary: factors of bandwidth u in ANY orthogonal grading => p >= sum_{f(t)-f(t')>2u} dim dim.
     u-local factors (any q): p >= N^2 Pr[L_1-L_2 > 2u]; banded factors: p >= (N-2w)(N-2w-1)/2; Bell with dim span S = r:
     p >= N^2 Pr[L_1-L_2 > 2r]. S' at u=0 ~ Theorem G for S_k; S' = robust Theorem G ("approximately symmetric terms").
 (b) odd-d Weyl-diagonal models: p >= (N^2-N)/2 (SWAP permutes qudit Bell states), tight [E d=3,5].
 Task 2 reductions: R1 delta_gen >= delta_Bell; R2 linear Bell corner => R_H(N^eta) <= N^{2-c eta} AND prob-rank_{N^{-c eta}}(H)
     <= N^eta via Fourier-sparse approximants (translation symmetry); R4 rigidity upper bounds do NOT transfer back; column
     relaxation p_Bell >= Rcol(q); coset-tight columns <= q^3 m^2/N^2.
NOT PROVED: any delta <= (1-eps) eta for general K (beyond line; impossible for class-B rank methods); delta = o(eta) even in Bell.
NUMERIC: S' within factor ~3-5 of best symmetric construction (k=1600); exact sector construction (101,35)@16, (452,56)@32;
 AP: (112,4)@16 general; AP unreliable (non-monotone), symmetric AP useless.
VERDICT: every structured family (site-symmetric, local, banded, bosonic/App B, clique-block, Weyl-diagonal) has a quadratic or
 worse corner, now PROVED for all S_k-symmetric/local ones. A linear corner (step 2) needs non-symmetric, highly non-local
 (locality >~ k sqrt(delta) in every factorisation), large-norm, exponential-growth structure; in the Bell submodel it would BE
 a new low-rank non-rigidity theorem for Walsh-Hadamard. Excluding it needs a beyond-line bound (rank methods cannot) -- in
 Bell a near-constant-fraction rigidity-type bound at polynomial rank for Fourier-sparse approximants of IP, far beyond the
 state of the art (constant-fraction rigidity of H known only up to rank ~log^2 N, F_3). Corner route: OPEN in both directions;
 recommend no further search inside symmetric/local/bosonic families.
Files: notes.md, sprime_check.py, band_check.py, sprime_asym.py(.out), sprime_d3.py(.out), hook_check.py, hook_min.py,
 sym_full_check.py, weyl_check.py, weyl_opt.py, sector_constr.py, ap_sym.py (+ ap_*.out), rcol.py, timing.py.
[N] N=16 general q=8: p=100 not converged (res .26); run killed at 10:18 (AP too weak to be informative).
