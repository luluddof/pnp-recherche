# G25-B2 — proof techniques beyond the line (brainstorm notes)
Tags: [P] short proof here (unreviewed), [E] exact check run here, [K] killed by argument, [H] heuristic.
Setting: F = swap on C^N (x) C^N, F = R + C, rank_Fl1 R <= p, osr C = rank_Fl2 C <= q. Goal: p q >= N^{2+eps} somewhere.

## 0. Two strategic remarks (read first)
(S1) WHERE to aim. HWY Prop 3.2 / 5.4 (G22a sect. 7): P(N) = min(p+q) >= N^{1+eps} gives exponential non-commutative circuit
  lower bounds for the permanent. So a bound beyond the line AT THE BALANCED POINT is a known-hard (major open) statement.
  An UNBALANCED bound (q = N^beta, beta small, p >= N^{2-beta+eps}) implies nothing about P(N) (feasibility is monotone, the
  balanced pair can still be cheap) and G18g says it already gives rho > 1/2.  => aim every technique at the extreme
  q = N^beta, beta -> 0 (target slope: T(theta) above 2 - theta near theta = 2; monomial model has 2 theta + T >= 4 there).
(S2) Rigidity dictionary. Theorem C = "untouched minor" argument (compress to a sub-block where the perturbation is cheap);
  the line pq ~ N^2 = the n^2/r barrier of untouched-minor rigidity bounds; Theorem B = the corresponding barrier.  Known ways
  past n^2/r in rigidity (Friedman's log, Lokam's bounded-coefficient / spectral rigidity) translate to: (a) adaptive choice
  of the compression (not covered by Theorem B, which is about one FIXED Phi), (b) NORM-controlled pieces (Lokam). (b) gives
  the first idea below and turns out to be EXACTLY the mechanism of the monomial bound.

## 1. Trace duality + "spectral spread" (generalises G22b Thm 3 to ALL decompositions)  [P + E]
 (a) For EVERY decomposition: N^2 = <F,F> = <F,R> + <F,C> <= ||R||_tr + ||Fl2 C||_tr <= p ||R||_op + q ||Fl2 C||_op
     (||F||_op = ||Fl2 F||_op = 1, Hoelder).  Also N^2 <= sqrt(p)||R||_F + sqrt(q)||C||_F.
 (b) In the monomial ST model, <F,R> = ||R||_tr exactly, ||R||_op = max T-class size <= g*(q) ~ sqrt(2q), same for C; so (a)
     IS G22b's Theorem 3(c) (N^2 <= sigma g*(tau) + tau g*(sigma)).  The intercalate key lemma is only an OPERATOR-NORM bound.
     [E] tdual.py on the g22b arrays: N = 8 Q8 array (8,8): <F,R> = <F,C> = 32, p op + q op = 64 = N^2 and the Frobenius form =
     64.000: BOTH inequalities are tight on the best monomial array.  (6;7,7): 42 vs 36; (7;8,8): 64 vs 49.
 (c) Bounded-spread model theorem (one line): if ||R||_op <= L and ||Fl2 C||_op <= L then p + q >= N^2/L.  Beats the line at
     the balanced point iff L <= N^{1-eps}.  Consequence: any decomposition with p = q = N^{1+delta} needs spectral norms
     >= N^{1-delta}/2; Appendix B (delta = 0.414) needs op norms >= N^{0.586}/2 (huge cancellation is forced).
 (d) SPREAD CONJECTURE: every feasible (p,q) has SOME decomposition (allowing p, q -> p N^{o(1)}, q N^{o(1)}) with
     ||R||_op <= N^{o(1)} (sqrt q + 1) and ||Fl2 C||_op <= N^{o(1)} (sqrt p + 1).  It would give the monomial curve
     2 theta + T >= 4 (balanced 4/3) for the GENERAL problem -- beyond the line, consistent with the circle (Appendix B
     satisfies the conclusion: 2 N^{sqrt2}(N^{0.707}) >> N^2).  Cancellation (border points, sigma_p n tau_q additions) is
     what can break it; the conjecture only asks for one good representative per orbit of the stabiliser GL x GL.
 Why beyond the barrier: Theorem B is norm-blind (rank only, same answer for F and for generic X); this is a Hoelder argument
 that uses ||F||_op = ||Fl2 F||_op = 1 with ||F||_F^2 = N^2 (F is unitary in BOTH flattenings: maximal trace norm / op norm
 ratio) and the decomposition-dependent spread L.  It already reproduces the 4/3 monomial bound, which is beyond the line.  Likely failure: optimal decompositions genuinely need op norms ~ N^{1-o(1)} (cancellation
 is intrinsic).  Already tried: no (G22d/G23d use distance, not the per-piece spread).  Test < 1 h: take the exact Clifford /
 Hurwitz / Hrubes decompositions (N = 2^k, G22c/G24a code) and the numerical (6,6)@6 of G24c; minimise max(||R||_op/sqrt q,
 ||Fl2 C||_op/sqrt p) over the stabiliser (a (x) b) . (b^-1 (x) a^-1) and over adding Z in sigma n tau (gradient descent);
 plot log L vs log N.  Prob. conjecture true: 25%; provable: 5%.

## 2. Disjoint-support (block) model: extend the 4/3 bound beyond monomial  [H + partial P]
 Model DS: R = sum_j a_j b_j^T with supp a_j pairwise disjoint and supp b_j pairwise disjoint (arbitrary complex entries);
 Fl2 C = sum_l vec(G_l) vec(G'_l)^T with disjoint supports.  Contains monomial ST (and non-symmetric, non-+-1 versions).
 Then every 4-tuple sees at most one R-term and one C-term (a "monomial-like" identity with real weights), and by (1a) the bound
 reduces to an op-norm lemma: an R-block with large overlap <a_j, F b_j> forces many C-blocks.  Observed mechanism [P sketch]:
 the off-diagonal nonzeros of Fl2(a_j b_j^T) = A_j (x) B_j sharing a row (x,u) must lie in ONE C-block (rows determine l), so each
 connected component of the tensor-product pattern supp A_j x supp B_j lies in one rank-one C-block; a dense k x k component
 then needs (swap_k) - A (x) B of rank one, impossible for k >= 3 by G19b Theorem 3 (k = 2 is the (1,1)@2 decomposition);
 sparse (matching-like) components give the intercalate count.  So: key lemma = intercalate counting + Theorem 3 on blocks.
 Why beyond line: same as monomial (counting, not rank).  Failure: medium-density components (between matching and dense) may
 escape both arguments; value limited (still a restricted model).  Already tried: only the strict monomial case (G22b, G23b).
 Test < 1 h: brute-force/SAT the DS support patterns at N = 4, 5 with LM on the values; check whether DS beats ST anywhere
 (G22b: exact (7,4)@5, (6,6)@5 use rank-2 R-blocks -> a "block-rank-r" DS variant is the right generality).  Prob 30%.

## 3. Friedrichs-angle form of the distance invariant (SDP-like)  [P]
 sigma_p + tau_q = U_{P,P'} (L_P + M_P'), L_P = {A : A = A P} (P rank-p projection), M_P' = {B : Fl2 B = P' Fl2 B}: LINEAR
 subspaces once (P, P') are fixed.  So c(N;p,q) = sup_{P,P'} ||Pi_{L_P + M_P'} F||^2, with ||Pi_L F||^2 = p, ||Pi_M F||^2 = q.
 Standard angle inequality: ||Pi_{L+M} f||^2 <= (||Pi_L f|| + ||Pi_M f||)^2 / (1 - gamma(P,P')) on the part orthogonal to
 L n M (plus the L n M part), gamma = cos of the Friedrichs angle.  Feasibility (c = N^2) forces 1 - gamma <= (sqrt p + sqrt q)^2/N^2:
 near-alignment of an "input-support" subspace with a "realigned-row" subspace, PLUS F correlating with L n M.
 Why beyond line: fixed (P,P') is exact linear algebra on F (not a flattening of X); the sup is over Grassmannians, not over a
 class-B map.  Failure: gamma -> 1 is easy (L n M contains |a(x)b><c(x)d|), so one needs the F-weighted angle; may just be a
 restatement.  Already tried: G23d Theorem A is the q = 1 case of this split (P only).  Test < 1 h: at N = 3, 4 take the LM
 optimisers of G22d, extract (P, P'), compute gamma and ||Pi_{LnM} F||^2; check whether c ~ (p+q)/(1-gamma) on all data.  Prob 8%.

## 4. Adaptive (Friedman-type) compressions: an unbounded factor beyond Theorem B  [H]
 Theorem B bounds one fixed Phi; choosing the compressing subspaces AFTER seeing (R, C) is outside it (as Friedman's untouched
 minor is adaptive).  Goal: p q >= N^2 log(N/q) (first unbounded-constant gain past the barrier; not an exponent).
 Failure (likely): for full-rank G_l, rank(pi_S G_l) = k for every S, and a generic p-dim range(R) has no special position; the
 combinatorial gain of rigidity comes from SPARSITY, which osr-pieces lack.  Works only if (1)/(2) first force structure.
 Not tried.  Test: at N = 4..6, for the numerical optimal decompositions, search over subspace pairs (S, S') maximising
 kN - rank Phi_S(R) - rank Phi_S(C) beyond the fixed-k value.  Prob 7%.

## 5. Positive (CP-sign) model  [P reductions]
 Over R, after the free sigma-symmetrisation (sect. 9), C = sum_l eta_l (G_l (x) G_l^T) and Pi*Fl1(F) = I.  The model eta_l = +1
 (C a completely-positive-type map, Choi PSD) is a GENUINE restriction: F(G (x) G^T) = M^T F M (M = 1 (x) G^T) is congruent to F,
 inertia (N(N+1)/2, N(N-1)/2), and the sum of q such forms must equal I + rank p (at most p negative eigenvalues).  Inertia +
 Weyl interlacing then constrain (p,q) without rank methods.  Failure: sums of congruent indefinite forms can be PD for q = 2
 generically; whether the Kronecker shape 1 (x) G^T prevents it is unclear; scale-free arguments impossible without norms.
 Test: numerically maximise lambda_min(sum_{l<=q} F(G_l (x) G_l^T)) / ||.|| at N = 3..6 to see the least q giving near-PD.  Prob 10%.

## 6. q-term spectral method (generalise G19b Thm 3)  [P partial]
 FK = I + E (rank p), K = sum A_l (x) B_l.  WLOG K commutes with F (cost 2x): K|S^2 = I + rank p, K|Lambda^2 = -I + rank p.
 If the A_l, B_l commute (simultaneously triangularisable), the eigenvalue on |ij>+|ji> and |ij>-|ji> coincide -> p >= N(N-1)/2.
 So feasibility needs the q-dim spaces A, B to be far from commutative; quantify by the dimension growth of words
 dim(AB)^k.  (FK)^2 = K'K in (BA)^k (x) (AB)^k.  Failure: generic q = 2 already generates M_N; no handle found.  Prob 5%.

## 7. Linear rank inequalities / network coding  [H]
 The crossing picture: x->v, y->u demanded; R = bottleneck {x,y}|{u,v} (log p), C = parallel channels x-u, y-v (log q).  The line
 theta + T >= 2 looks like a cut-set bound; for LINEAR codes, Ingleton / DFZ inequalities beat cut-set bounds.  Build a subspace
 configuration from (range R, row R, span G_l, span G'_l, ...) and apply Ingleton.  Failure: no natural configuration with
 F-specific ranks found in 10 min.  Prob 3%, value high, cost high.

## 8. Higher-degree covariants (continue G21d)  [—]
 Degree-8 separating weights at N = 3 for J(2,2).  Value low (one witness at N = 3 says nothing about exponents).  Prob of a
 scalable template 2%.

## 9. KILLED in this session
 * "G' = G^T symmetric pieces" as a wider model [K, P]: the involution sigma: X(x,y,u,v) -> X(v,u,y,x) fixes F and preserves
   rank (Fl1 -> Pi Fl1^T Pi) and osr; averaging gives a sigma-invariant decomposition at (2p, 2q); sigma-invariant means Pi Fl1(R)
   is complex symmetric, Takagi gives R = sum h_j h_j^T-type (G'_l = G_l^T form).  Over C the symmetric model = general model up
   to a factor 2; over R it is the SIGNED model (sect. 5 is the unsigned restriction).
 * Asymptotic spectrum / quantum functionals [K]: feasibility implies the restriction F <= (GHZ_p (x) EPR_xy (x) EPR_uv) (+)
   (GHZ_q (x) EPR_xu (x) EPR_yv); every spectral point gives N^{sum theta} on every EPR pair product, hence only p + q >= 1.
 * Hopf-Stiefel / topology [K]: these bound nonsingular bilinear maps, which exist with n ~ 2N; capped at linear, useless.
 * Entanglement monotones [K as a separate route]: Schmidt RANK = flattening (class B); entropic monotones are not SLOCC-monotone
   (non-unitary G_l amplify small Schmidt coefficients); SLOCC invariants = covariants = sect. 8.  Mixed-input trick
   (rho of rank r: rho^{T_B} = osr p + rank qr) only weakens.
 * Twirl / Fourier / linear witnesses [K]: max <Y,F> under ||Y||_op, ||Fl2 Y||_op <= 1 is N^2 = trace duality (sect. 1); no new
   information without norm control.
 * Dimension count on rank-r vectors in ker R (rowspace(W) in G*colspace(W)) [K]: optimising r gives p >= N^2/(4q) = line again.

## TOP 3 (probability x value / cost)
 1. Spread conjecture test (sect. 1d): measure min op-norm spread of known non-monomial decompositions vs N.  Cheap, decides
    whether the monomial mechanism can be general.
 2. DS / block model (sect. 2): turn "intercalate counting + Theorem 3 on dense blocks" into a p^2 q >~ N^4/poly bound.
 3. Friedrichs-angle form (sect. 3) on the G22d optimisers, aimed at q = N^beta small.
Files: notes.md, tdual.py (trace-duality check on g22b arrays), sigma.py ([E] sigma fixes F, preserves both ranks, symmetrises).
 [E] sigma.py: sigma F = F; ranks (2,2), osr (3,3) preserved; Pi Fl1(R + sigma R) and Pi' Fl2(C + sigma C) symmetric.
