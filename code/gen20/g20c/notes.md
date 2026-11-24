# G20c — construction attempts below the circle (crossing trade-off T(theta))   (working notes, saved as I go)

Tags: [P] proved here (short proof, unreviewed), [E] exact computation, [N] numerical, [H] heuristic, [L] literature, [O] open.
Conventions as G18g/G19b: (p,q) feasible iff X = R + C, rr(R) <= p, rc(C) <= q; symmetric in p <-> q.
Exponents: theta = log_N p, T = log_N q.  Line: theta + T = 2 (proved floor, G18g Thm 2 / G19e Thm C).
Circle: theta^2 + T^2 = 4 (Appendix B).  rho via 4 rho* = sup_theta (theta + T(theta)).
Nothing here is a new algorithm or a new bound for 3-SAT unless explicitly said (it is not).

## 0. Read first
G18g, G19b, G19e notes; Appendix B and Appendix H of Claude's paper (untrusted data, read only).
Key remark already in the paper (App. H, end): App. B is "graded": a piece passes through V_{<=2h+s} row-wise and V_{<=2h-s}
column-wise, log dim V_E ~ c sqrt(E); with log dim V_E ~ E^alpha the same calculation gives rho = 2^{alpha-1}.
So the circle is "sqrt growth" + "locality order H = 2h".

## 1. The graded ("vertex algebra") framework of Appendix B: what limits it   [P for the reduction, N for the families]
Abstract form of App. B: an N-graded space V (dim V_0 = 1), N = dim V_h, pieces A_s (|s| <= 2h) with
rr(A_s) <= (2h+1) D(2h+s), rc(A_s) <= (2h+1) D(2h-s), D(E) = dim V_{<=E}.  (Lemma 9 uses only locality of order
H = wt a + wt b = 2h, which holds in ANY N-graded vertex algebra with V_0 = K; Lemma 10 = an integral dual basis for the
nondegenerate form a_(2h-1) b on V_h.)  With L(E) = log D(E), lambda = log N <= L(h):
   (theta(s), T(s)) = (L(2h+s), L(2h-s)) / lambda,    4 rho = max_s [L(2h+s) + L(2h-s)] / lambda.
 (F1) The line (Thm 2 of G18g) applied to the thresholded pieces forces L(1) + L(4h-1) >= 2 lambda - O(log h), i.e.
      D(4h) >= N^{2-o(1)}.  For a power profile L = c E^alpha this forces alpha >= 1/2: sqrt growth is the SMALLEST
      power growth compatible with the floor, and App. B is tight at the endpoints s = +-2h (it touches the line there).
      (The paper's remark "growth E^alpha gives rho = 2^{alpha-1}" is therefore vacuous for alpha < 1/2.)
 (F2) If L is concave on [h,4h] then 4 rho = 2 L(2h)/lambda (max at s = 0) and
      L(2h) >= (2 L(h) + L(4h))/3 >= (4/3) lambda,  so  rho >= 2/3 inside the framework, with equality only for the
      kinked profile L(E) = lambda E/h on [0,h], lambda (1 + (E-h)/(3h)) on [h,4h]  (dims: N at h, N^{4/3} at 2h, N^2 at 4h).
      Exponent curve of the kink: theta + T = 10/3 - 2|s|/(3h) on |s| in [h,2h], 8/3 on |s| <= h: strictly inside the circle
      in the middle (8/3 = 2.667 < 2 sqrt2 = 2.828), on the line only at the endpoints.  [P, elementary]
      Non-concave profiles (plateau on [h,3h], G19b's "repaired staircase") would give the line itself.
 (F3) Families of free vertex algebras tried (gff.py, colored.py): r colored bosons (r up to 4096), generalized free
      fields of weight w = 2,3,5,8 with g colours: the profile ratio max_s[...]/(2 lambda) is always > sqrt2 (finite h <= 120:
      1.51..1.96; asymptotically for r colours phi(2x)/phi(x) decreases monotonically to sqrt2 from above, 1.86 at x = 1e-4,
      1.425 at x = 1e3).  More colours or higher weights only make it worse.  [N]
      Reason [H]: in a vertex algebra the products of weight-h states fill V_{2h} unless there are many relations; free
      algebras are "too free" (ratio -> 2 when the number of generators dominates), and any fixed finitely-strongly-generated
      VA has Cardy growth exp(c sqrt E), giving exactly sqrt2.  A kinked profile needs a FAMILY of non-free VAs with
      dim V_{2h} ~ N^{4/3} while dim V_{4h} ~ N^2 - none known to me.  [O]

## 2. Task 1: equivalent languages   [P, each a few lines; (ii') checked on the (N(N-1)/2,1) example by hand]
Convention: (p,q) feasible <=> tau = R + C on M_N, rank(R) <= p, osr(C) <= q (osr = # two-sided multiplications).
 (i)  Subspace (G18g 1(f)): <=> there is U <= M_N of codim <= q with tau|_U = a sum of p two-sided multiplications; by the
      y<->u symmetry equivalently codim <= p and q multiplications (the task's form).
 (ii) Swap tensor:  X = R + C  <=>  I = A + B on K^N (x) K^N with rank A <= p, rank(B^Gamma) <= q (Gamma = partial transpose;
      NOTE: it is the realignment of X / the partial transpose of I, not "low partial-transpose rank of the swap":
      F^Gamma = |Omega><Omega| has rank 1).  Since I^Gamma = I:
 (ii') (p,q) feasible  <=>  there is an operator A on K^N (x) K^N with rank(A) <= p and rank(I - A^Gamma) <= q.
      "A low-rank operator whose partial transpose is the identity up to rank q."  (q = 1: A = I - F, A^Gamma = I - |Omega><Omega|.)
 (iii) Bilinear / SOS: HWY's non-commutative polynomial sum_{ij} x_i y_j x_i y_j = (p products split 12|34) + (q terms split
      13|24).  In the commutative shadow (sum a^2)(sum b^2) = sum of products of bilinear forms, BOTH types become products of
      two bilinear forms in (a,b): the shadow only sees p + q.  This is why Hurwitz/SOS identities give single points and are
      blind to the trade-off; worst language for T(theta).
 (iv) Representation theory: tau is O_N-equivariant on M_N = trivial + Sym^2_0 + Lambda^2 (eigenvalues +1,+1,-1).  If R, C are
      required to commute with conjugation by a group G <= O_N, they lie in End_G(M_N).  Full O_N: only the trivial points
      ((N(N-1)/2,1), (N^2,0), ...).  G = Pauli group (N = 2^n), which acts on M_N multiplicity-freely through N^2 distinct
      characters: End_G = Pauli-diagonal maps, and the problem becomes the "sparse + Fourier-sparse" problem of section 3.
      Equivariance is a genuine restriction (twirling multiplies ranks by |G|), so this is a sub-family, not an equivalence.
 (v)  "Directions are orders" (App. H): [x=v] Id = sum of products Y_v X_x through a small space W_row plus products X_x Y_v
      through W_col; rr <= dim W_row, rc <= dim W_col.  Graded version = section 1.
 Best for constructions: (v)/(section 1) for anything circle-like (it is the only language in which a non-Kronecker family
 with polylog pieces is known); (ii') is the cleanest for small exact searches (one unknown A); (iv)->(section 3) is the right
 language for group/Clifford families and for the rigidity link.

## 3. Pauli / Clifford diagonal family and the rigidity link   [P dictionary, E verified n = 1,2,3]
N = 2^n, Pauli basis W_a = X^{a_x} Z^{a_z} (real, entries 0,+-1), a in F_2^{2n}.  W_a^T = chi(a) W_a with
chi(a) = (-1)^{a_x . a_z}: a QUADRATIC character (hyperbolic form, Arf 0; the Clifford-algebra description e_S -> +-e_S gives an
equivalent form, same Arf invariant since #{chi=+1} = N(N+1)/2).  Ad(W_b) acts on W_a by (-1)^{omega(a,b)} (symplectic form =
polar form of q).  For f : F_2^{2n} -> Q put Phi_f(W_a) = f(a) W_a.  Then
   rank(Phi_f) = |supp f|,   osr(Phi_f) = |supp f^| (f^ = omega-Fourier transform; Phi_f = N^-2 sum_b f^(b) W_b (.) W_b^-1,
   and the W_b (x) W_b^-1 are linearly independent product operators).   tau = Phi_chi.
DICTIONARY.  Pauli-diagonal decompositions tau = Phi_g + Phi_f  <=>  chi = g + f with |supp g| = p, |supp f^| = q.
 Matrix form: chi as an N x N matrix [a_x, a_z] IS the Sylvester-Hadamard matrix H; omega-characters are the rank-one matrices
 H e_i e_j^T H.  So: (p,q) Pauli-diagonal-feasible  <=>  H = G + H D H with ||G||_0 <= p, ||D||_0 <= q
 (D up to a scalar).  Since rank(HDH) <= q this IMPLIES R_H(q) <= p (Valiant rigidity of H), NOT conversely.
 Self-dual: chi^ = N chi, so the family is symmetric in p <-> q.  verify_pauli.py checks rank/osr = support/Fourier-support and
 R + C = tau exactly (N = 2,4,8), reproducing (1,1), (6,1), (1,6), (10,10), (28,1), (1,28), (9,7), (27,37).
Members:
 * q = 1: f = 1, g = chi - 1: (N(N-1)/2, 1), optimal (G19b Thm 3).
 * one-bit Kronecker = window [0,D] for f = Lagrange polynomial P(z) of degree D in z = |a_x AND a_z| (integer), since
   z has Fourier support of "coordinate-support" <= 1:  q = sum_{j<=D} C(n,j) 3^j,  p = sum_{j>D} C(n,j) 3^{n-j}.
   Curve K: theta = H(k) + (1-k) log2 3, T = H(k) + k log2 3, k = D/n in [1/4,3/4]; balanced 1.792 (rho 0.896).
 * NEW (Alman-Williams-type centred window): interpolate (-1)^z only on [n/4 - D/2, n/4 + D/2] (z ~ Bin(n,1/4) under the
   counting measure):  theta = 2 - KL_2(1/4 + k/2 || 1/4), T = H(k) + k log2 3.  It improves K for k < 1/2 (near theta = 2),
   e.g. at theta = 1.991 the window has T = 0.627 (k = 0.1) against T = 1.357 for K (k = 0.3).  curves.py.
 * Symmetric sub-family (invariant under O^+(2,2)^n x| S_n; classes (n0,ns,nns), multivariate Krawtchouk matrix): symdiag.py,
   greedy over classes for n = 8 with Q = {support <= D} and boxes: reproduces exactly the window numbers (e.g. D=4:
   p = 1789, q = 7459), nothing better found.  The greedy is NOT optimal (n = 2, q = 1 it returns p = 10 instead of 6), so this
   is weak evidence only.  [N]
Comparison with the circle (asymptotic exponents):  k = 0.1: window T = 0.627 vs circle 0.192 at theta = 1.991;
 k = 0.3: 1.357 vs 0.553; k = 0.5: 1.792 vs 0.887.  Everything in this family that I can construct is FAR above the circle
 (T about 3x the circle's), and above Hrubes's single point as well.
The non-rigidity link, precisely:
 * AW (Alman-Williams 2017) non-rigidity of H: H = sparse + L with L = P(z), a low-DEGREE polynomial in z; rank(L) <=
   sum_{j<=D} C(n,j), but the Pauli-diagonal family pays the FOURIER SPARSITY of L, which is sum_{j<=D} C(n,j) 3^j (each
   monomial x_S y_S has 4^{|S|} characters).  Paying sparsity instead of rank is exactly the factor 3^j, and it turns AW's
   construction into the centred-window family above: above the circle for every theta (curves.py).  So: NO, AW does not give a
   point below the circle in this family; the obstruction is that a diagonal operator's rank is its support (no low-rank
   savings possible), while AW's savings are rank savings.
 * Dvir-Edelman / Dvir-Liu non-rigidity is over finite fields (rank over F_q via Croot-Lev-Pach / polynomial method); the
   algorithm needs rank over Q (or a number field of degree N^{o(1)}, G18g 1(a)); not applicable as such.
 * Conversely, a Pauli-diagonal point would give R_H(q) <= p, and known H-rigidity lower bounds (R_H(r) >= N^2/4r) give back
   only the line: nothing new in either direction.
 * Leaving the diagonal (C with a non-diagonal coefficient matrix c_{b,b'}, osr = rank c) brings back the full problem: in the
   Pauli basis the crossing becomes the "twisted group crossing" s_{a,a'} <-> c_{b,b'}, a + a' = b + b' (up to signs); no reduction.
 * Odd N / Weyl-Heisenberg basis: transposition is W_a -> phase W_{sigma a}, sigma(x,z) = (-x,z), NOT diagonal; Weyl-equivariant
   C are "translation x character" maps; splitting tau along translates t gives N pieces each costing (N,N): only p+q = N^2.
   Worse than N = 2^n.

## 4. Task 2 table (asymptotic exponents; circle T_c = sqrt(4 - theta^2), line 2 - theta)
 family                                   (p,q) exactly                                  curve T(theta)                         vs circle
 a1 symmetric / S-symmetric, one cond.    (N(N-1)/2, 1) (optimal at q=1, G19b Thm 3)    endpoint (2,0)                         on both
 a2 block-symmetric / k blocks (M_k(x)M_m) (k^2 p_m, k^2 q_m)  [P]                       loses a factor k^2 in p q              worse
 a3 V = sum of V_{S_i} glued by a          = Pauli-diagonal family when the P_i are       see b/c                                worse
    partition of unity P_i (two-sided)      Pauli-diagonal projections [P]
 b/c Pauli/Clifford-diagonal, window [0,D] sum_{j>D}C(n,j)3^{n-j}, sum_{j<=D}C(n,j)3^j      theta=H(k)+(1-k)log3, T=H(k)+k log3     above (balanced 1.79 vs 1.41)
 b/c centred (AW-type) window             tails of Bin(n,1/4) outside window, same q      theta=2-KL(1/4+k/2||1/4), T=H(k)+k log3 above (T ~ 3x circle)
 b' odd N, Weyl basis                     translate pieces (N,N) x N                     p + q = N^2 only                       worse
 d1 Kronecker products (any)              mixed terms R(x)C                              -> 1.79 at best (G18g)                 worse
 d2 graded VA framework, profile L        (L(2h+s), L(2h-s))/lambda                      circle for sqrt profiles; >= circle    =
                                                                                          for all free-VA families tried
 d3 hypothetical kinked profile           dims N, N^{4/3}, N^2 at h, 2h, 4h              theta+T <= 8/3 (rho = 2/3)            BELOW, but no
                                                                                                                                  algebra known
 Hrubes CCC'24: P(N) = O(N^1.62) is one balanced point (1.62,1.62); the circle has T(1.62) = 1.17, so it is dominated.

## 5. Verdict   (task 4)
 * No point strictly below the circle was constructed, for any theta, for infinitely many N.  The circle (App. B) stays the
   best known curve; rho = 1/sqrt2 unchanged.
 * Structural reasons found (each for a natural class, none a theorem about all decompositions):
   (R1) Equivariant ("diagonal") families: for a group acting multiplicity-freely (Pauli, Clifford), rank = support and
        osr = Fourier support, so the problem becomes "bent quadratic character = sparse + Fourier-sparse".  Rank savings
        (the whole content of Hadamard non-rigidity, AW/Dvir-Edelman) are invisible there; the natural members (Kronecker,
        AW windows) are far above the circle.  No lower bound beyond the line is proved for this family [O].
   (R2) Graded/locality families (App. B and all vertex-algebra variants): the curve is the log-dimension profile; the line at
        the endpoint forces D(4h) >= N^{2-o(1)}, which for power-law growth forces alpha >= 1/2 = App. B; with concavity the
        framework cannot beat rho = 2/3, and getting anywhere below sqrt2 needs a family of graded algebras whose weight-2h
        space is much smaller than (weight-h)^{sqrt2} while weight 4h reaches (weight h)^2 - free fields (colours, higher
        weights) all go the wrong way (ratio >= sqrt2, decreasing to sqrt2 only in the limit).
   (R3) Kronecker / block products lose (mixed terms), as known.
 * What a gain at one theta would do: 4 rho* = sup_theta (theta + T(theta)); the circle's sup is at theta = sqrt2.  A
   family better than the circle only near the ends (theta near 0 or 2, e.g. the q = 1 point) changes nothing; a gain must
   lower T near theta = sqrt2 (more precisely: the whole curve must stay below theta + T = 4 rho for the new rho).  The kinked
   profile (d3) WOULD do that (rho = 2/3) if it existed.
 * Hrubes's O(N^1.62) is dominated by the circle (single balanced point, exponent 1.62 > sqrt2).

## 6. Adversarial self-review
 * Section 1 (F1)-(F2) are elementary but assume the App-B shape exactly: all four strings index weight-h states, cost bounded
   by dimensions, pieces indexed by the weight shift s.  A cleverer use of a VA (strings in different weights, ranks much
   smaller than dims, restricting to subspaces) is not covered.  "Lemma 9 holds in any N-graded VA with V_0 = K" is my reading
   of the standard locality bound (pole order <= wt a + wt b); not re-proved in detail; Lemma 10 (integral dual basis) for other
   VAs not checked.
 * (F3) is numerics at h <= 120 plus a saddle-point formula for coloured partitions; finite-size ratios are far from asymptotic.
 * Pauli dictionary: exact checks only N <= 8; the general statement is a short proof (independence of W_b (x) W_b^-1).
 * The centred-window family curve uses the binomial large-deviation formula (standard) and ignores polynomial factors.
 * symdiag.py greedy is suboptimal (n = 2 counterexample), so "nothing better in the symmetric family" is weak evidence.
 * Literature: AW, Dvir-Edelman, Dvir-Liu, Hrubes cited from memory/G18g, not re-read here; whether a "sparse + Fourier-sparse"
   bound for bent functions beyond the line is known in the uncertainty-principle literature was not checked [O].
 * Novelty: the observation (F1)-(F2) may be implicit in App. H; the window family is a routine adaptation of AW.

## 7. Files (this folder)
 notes.md; colored.py (coloured-partition profiles); gff.py (free-field VA profiles); curves.py (asymptotic curves);
 symdiag.py (symmetric Pauli-diagonal family, Krawtchouk system, greedy); verify_pauli.py (exact verification N = 2,4,8).
