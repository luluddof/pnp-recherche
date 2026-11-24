# G21a — can a module / twisted module / degenerate ground space lower rho below 1/sqrt2?   (working notes)
Tags: [P] proved here (short, unreviewed), [E] exact computation, [N] numerical, [H] heuristic, [L] literature from memory, [O] open.
VERDICT: NO construction. No candidate is below the circle, already at the level of exact dimension profiles.
Nothing here is a new algorithm or a new bound for 3-SAT.
What rho would give (Theorem A of the paper, alpha' = log2(1/rho); Def. 7 needs 1/2 <= rho < 1):
  rho = 1/sqrt2: 2^{O(km/(log m)^{1/2})} (current);  rho = 2/3: alpha' = 0.585, i.e. the POWER of log m rises from 0.5 to 0.585
  (not just a constant);  rho = 1/2: alpha' = 1, 2^{O(km/log m)}; in n for fixed k: 2^{O(n loglog n/(log n)^{alpha'})} (Thm A(c)).

## 1. Abstract theorem extracted from Appendix B   [P modulo standard vertex-algebra facts; E for Fock modules]
Data.
 (A1) V a vertex (super)algebra over Q, M a V-module (or a g-twisted module, g of finite order T), with ANY grading
      M = (+)_{E >= 0} M_E (E in (1/T')Z, dim M_E finite) for which homogeneous modes shift degree: a[r] M_E c M_{E+r}.
      (The grading need not be the conformal one; only locality orders are intrinsic.)
 (A2) states a_x, b_v in V (x, v in [N]) and an INTEGER H with (z-w)^H [Y_M(a_x,z), Y_M(b_v,w)]_{+-} = 0 for all x, v
      (true as soon as a_x (n) b_v = 0 in V for n >= H; for twisted modules the exponent of (z-w) stays the same integer,
      only the powers of z, w become fractional [L: twisted Jacobi identity, not re-derived]).
 (A3) labels: N basis vectors y of M_{E0} with  <u| o(a_x (H-1) b_v) |y> = d [x=v][y=u], d != 0
      (o = weight-0 mode; sufficient: a_x (H-1) b_v = d [x=v] 1.  In App. B: Lemma 9(c) + Lemma 10, H = 2h, E0 = h, M = V.)
 (A4) integrality/uniformity: a Z-form (integer mode matrices of poly(lambda) bits, d integer, tables in 2^{O(lambda)}) (O1, Lemma 10).
Conclusion.  With B = 2 E0 + H and Q_s = sum_{i=0}^{H} (-1)^i C(H,i) a[..+i] b[..-i] (weight-0 coefficients of
 (w1-w2)^H Y(a,w1)Y(b,w2)), the pieces A_s = <u|Q_s(a_x,b_v)|y> satisfy sum_s A_s = d X, there are <= T'(B+1) of them, and the
 piece with "row weight" x factors row-wise through (+)_{E in [x-H,x]} M_E and column-wise through (+)_{E in [B-x-H, B-x]} M_E:
      (theta, T)(x) = (L_M(x), L_M(B-x))/lambda,   4 rho lambda = max_x [L_M(x) + L_M(B-x)] + O(log),   L_M = log dim M_{<=E}.
 N is bounded by BOTH: N <= dim M_{E0} (labels) and N <= N_V(H) := max number of dual pairs of locality order <= H;
 for V of CFT type (V_0 = Q1, self-dual) N_V(H) <= dim V_{<= H/2} in the conformal (charge-symmetric) grading  [(V')].
Proof = proof of Lemma 9 read on M (weak commutativity + weak associativity on modules); budget: the a-first products pass
 through weights <= E0 + s + h_a, the b-first ones through <= E0 + H - h_a - s; the sum is 2E0 + H whatever the grading.
What does NOT fit [P: the argument breaks]: generalized vertex algebras / intertwining operators / parafermions, where the
 exponent of (z-w) is fractional: the binomial series of (w1-w2)^H is infinite and the a-first window is unbounded.
 Non-vacuum normalisation c = a (H-1) b with o(c) scalar on the labels shortens H by wt c only (c = omega: 2), as G18h said.
Exact check [E] (appb_mod.py, appb_mod.out): Heisenberg Fock modules pi_gamma (zero mode gamma = 0,1,2,3), labels at weight
 E0 in {2,..,6} != h, h in {2,3,4}: sum_s A_s = h![x=v][y=u] EXACTLY in all 7 runs, and every exact threshold rank is <= the
 window dimension P(E0+h+s0-1), P(E0+h-s0) (equal up to a few units), i.e. budget B = 2E0 + 2h as stated.
 (Mixed-weight labels in M_{<=E0} should work the same way [P sketch, not run].)

## 2. Profile optimisation, corrected   [P, all in the dimension accounting]
Constraints: L_M(E0) >= lambda;  rank count / line theorem: L_M(x) + L_M(B-x) >= 2 lambda - O(log) for every x, so L_M(B) >= 2 lambda.
 (G) Ground-space penalty:  4 rho lambda >= L_M(B) + L_M(0) >= 2 lambda + L_M(0),  i.e.  rho >= 1/2 + L_M(0)/(4 lambda).
     A degenerate ground space is a pure LOSS in the accounting; labels in the ground space (E0 = 0, L_M(0) >= lambda) give
     rho >= 3/4.  (An inert tensor factor would not be charged by true ranks, but it is also useless.)
 (C) Concave L_M:  rho >= 1/2 + H/(4(E0+H))   (interpolate between (E0, lambda) and (B, 2 lambda) at B/2).
     H = 2E0 (App. B shape: M = V, E0 = h, H = 2h): 2/3.   E0 = 0: 3/4.   H/E0 -> 0: 1/2.
 Where 2/3 uses concavity: only at the three points h, 2h, 4h, AND it uses M = V, E0 = h, H = 2h.  Without concavity
 (monotone only): 1/2 (plateau = repaired staircase).  Without tying H to E0: -> 1/2 even for concave L (G18h's scenario (1)).
 Coordinator's lead, checked: the affine profile L = (2/3 + E/(3h)) lambda is AN extremal for 2/3 (cost 8/3 at every x), but
 so is the kink with L(0) = 0; the intercept 2 lambda/3 is the largest ground space that (G) ALLOWS at rho = 2/3, not something
 the bound needs.  What is needed is the slope lambda/(3h) on [h,4h].  "Degenerate ground space helps" is false as stated.
 (L) Local-exponent lemma [P]: for S(E) = min_b [log Z(b) + bE],  d log S/d log E < 1/2  <=>  d/db [b log Z(b)] > 0.
     For M = V a gain needs exponent < 1/2 somewhere on [h,2h] and > 1/2 on [2h,4h] (the line forces L(4h) >= 2L(h)).

## 3. Candidates: exact profiles (profiles.py/.out, bigh.py/.out, locexp.py/.out) and exact small cases
ratio := max_x [L_M(x)+L_M(B-x)]/(2 lambda)  (circle sqrt2 = 1.4142, rho = ratio/2), N = dim V_{<=h}, E0 least with dim M_{<=E0} >= N.
 object (V ; M)                                      profile                                   best exact ratio (h <= 30)        small-case exact check
 Heisenberg rank 1 ; V (App. B)                      p(E)                                      1.479 (h=119)                      identity + ranks, h<=7 (G20b)
 Heisenberg ; Fock module pi_gamma, labels at E0     same                                      same as App. B                     [E] appb_mod.py: identity exact, ranks = windows
 Heisenberg rank d ; Z2-twisted (half-odd modes)     prod (1-q^{k-1/2})^{-d}                   1.61 (d=2) .. 1.88 (d=64)          not run
 lattice-type twisted, ground 2^{d/2}                2^{d/2} x the same                        1.62 .. 2.04                       not run
 n fermions NS ; NS                                  prod (1+q^{k-1/2})^n                      1.4395 (n=32), -> sqrt2 from above  not run
 n fermions NS ; Ramond (ground 2^{n/2})             2^{n/2} prod (1+q^k)^n                    1.4427 (n=32), 1.4279 (n=16,h=160)  [E] n=8, H=1: = 8-square identity
 complex fermions ; regradings bc(1/4,3/4), bc(0,1)  spectral flows                            1.443 .. 1.461                      not run
 symplectic bosons (c=-n) ; V                        prod (1-q^{k-1/2})^{-2n}                  1.476 (n=8) .. 1.537 (n=64)         not run
 symplectic fermions (c=-2n) ; V                     prod (1+q^k)^{2n}                         1.501 .. 1.583                      not run
 "finite filling" exterior profile log C(2n,2E)      not a module of a vertex algebra          1.52 at the largest h the line allows   -
 In every row the line check min_x cost/(2 lambda) >= 1.03 holds (consistent with the Lean floor).  NO row is below sqrt2.
Ramond, smallest instance [E, ramond8.py]: n = 8, a_x = b_x = psi^x (H = 1), labels = chiral ground space (dim 8):
 {psi^x_0, psi^v_0} = [x=v]: R = <u|g_v g_x|y>, C = <u|g_x g_v|y>, R + C = 2X exactly, rr(R) = 8, rc(R) = 64, rr(C) = 64, rc(C) = 8.
 This IS the eight-square identity of the paper's Lemma 2 (App. H (74)): pair (8,8) on the line, but as a piece 9 of 12 bits,
 rho = 3/4 = exactly bound (G) with L_M(0) = lambda.  So the Ramond lead reproduces the old Hurwitz-Radon construction.
Why Ramond does not help [N + H]: with log Z = a/b - (c/24 - Delta_M) b + sigma (one-term modular expansion; sigma = 0 for
 fermions), S_V(h) = A sqrt(h - c/24), S_M(E) = A sqrt(E + Delta_M - c/24); the Ramond MODULE profile alone is the wanted
 "affine with intercept" shape (exponent < 1/2 at every energy: locexp.out, row "Ramond": g > 0 everywhere) but the number of
 states is counted in V: the midpoint of the budget sits at cylinder energy u_M + u_V + c/24 = 2u + c/24 for EVERY module, so
      ratio^2 = 2 + c/(24u) > 2   (u = (lambda/A)^2),  independent of M.
 bigh.out: exact ratios for n = 8, 16, h up to 160 decrease monotonically to sqrt2 from above (1.4313, 1.4273 at h = 160;
 leading prediction 1.4146, 1.4150; the gap is the E^{-1/4} prefactor).  Ramond is slightly WORSE than NS at each size.
General one-term ("Cardy") statement [P given the asymptotic form, H for actual dimensions]: with p = (lambda - sigma_V)/A,
 q = (lambda - sigma_M)/A, cost/2 = A sqrt(p^2 + q^2 + c/24) + sigma_M; if c >= 0 and sigma_V <= sigma_M <= 0 (all unitary
 rational V, all lattice and free-fermion algebras, any module) then ratio >= sqrt2.  A gain needs c < 0 with
 |sigma_V| in a narrow window (for M = V: 2p - sqrt(4p^2-3|w|) <= |sigma_V| < (sqrt2 p - sqrt(2p^2-|w|))/(sqrt2-1), w = c/24),
 whose limit is exactly the affine extremal rho = 2/3.  Symplectic fermions come nearest (|sigma_V| = n log 2 against
 A sqrt|w| = 0.74 n) but only where the expansion is invalid; exactly they are at 1.50-1.58.
Free fields, M = V [N on a grid of b, analytic at both ends]: every tower allowed in a vertex algebra (bosons of weight
 1/2, 1; fermions of weight 1/2, 1; higher weights are more negative) has d/db[b log Z] < 0 at ALL temperatures
 (locexp.out), so by (L) the local exponent is > 1/2 everywhere and L(2h) > sqrt2 L(h): this upgrades G20c (F3) from
 "tried" to a structural statement for free-field algebras.  Towers with g > 0 (fermions of weight < 1/2, Ramond zero
 modes) exist only in modules/regradings, where they come paired with a conjugate tower that cancels them, or via (V').

## 4. Obstruction, as precisely as I can state it
 (O1) [P, accounting] rho >= 1/2 + L_M(0)/(4 lambda); labels in a ground space: rho >= 3/4 (attained: Ramond n = 8 [E]).
 (O2) [P, accounting] for a module M != V the budget is B = 2E0 + H with lambda <= min(L_M(E0), L_V(H/2)+1): a module can
      only help if it grows SLOWER than V per unit of locality order (H << E0), the opposite of a large top space.
 (O3) [N/H] in every family computed, and in the one-term modular regime for c >= 0, sigma <= 0, the module-independent
      identity "midpoint energy = 2u + c/24" gives ratio >= sqrt2.  Not a theorem about all vertex algebras: non-unitary
      families with c < 0, c_eff > 0 and a wide one-term regime (large gap in the dual channel) are NOT excluded [O];
      the Yang-Lee series M(2,2k+3) fails because its dual-channel gap is ~1/k.
 Exact unproved steps: twisted locality [L]; N_V(H) <= dim V_{<=H/2} for non-CFT-type V; g < 0 for all b analytically;
 passage from Legendre entropy to actual dimensions in scaling families.

## 5. Adversarial self-review
 * No exact rank computation was done for fermionic or twisted pieces beyond n = 8, H = 1; rows marked "not run" are
   profile-level only.  Since none beats the circle at profile level and ranks <= dimensions, this cannot hide a gain
   UNLESS true ranks are far below window dimensions (G20b found them equal for App. B; not tested for modules with ground space).
 * Ratios are finite-size (h <= 30, lambda <= 170) plus n = 8, 16 to h = 160; the claim "-> sqrt2 from above" rests on the
   one-term formula, checked only in trend.
 * The grading freedom in (A1) and (V') are my reading; "N_V(H) is grading-free" was argued, not machine-checked.
 * E8/Leech-type holomorphic lattice algebras were not computed (covered only by the heuristic c > 0, sigma = 0 statement).
 * The bound (C) -> 1/2 as H/E0 -> 0 is an accounting possibility, not a construction.

## 6. Files (this folder)
 notes.md; ramond_so8.py/.out (Ramond n=8, so(8) currents, H=2, exact); locexp2.out; profiles.py/.out (exact profiles); bigh.py/.out (large h, fermions); locexp.py/.out (local-exponent scan);
 appb_mod.py/.out (exact pieces on Heisenberg Fock modules, labels at any weight); ramond8.py/.out (Ramond n = 8 = 8 squares).

## 7. Reply to the reviewer G21e (relayed by the coordinator): do Lemma 9 / Lemma 10 extend to modules and direct sums with a
##    degenerate lowest-weight space, and do exact ranks follow a plateau profile?
 * Lemma 9 EXTENDS [P: it is locality + associativity, valid on every module, every direct sum, and (integer exponent) on
   twisted modules [L]].  Exact confirmations [E]: Heisenberg Fock modules pi_gamma (7 runs, sect. 1); Ramond n = 8 with
   H = 1 (ramond8.py) and with H = 2 (ramond_so8.py: a_x = b_x = psi^i psi^j, N = 28 so(8) currents, labels = the
   16-dimensional ground space on which the zero modes act as the Clifford algebra): sum_s A_s = -4 [x=v][y=u] EXACTLY.
 * Lemma 10 is a statement about V alone (a[-h] b in V_0); it does not see the module.  Nothing to extend.
 * So "unknown" is resolved: the lemmas hold.  What fails is the PROFILE: the pieces only visit the V-submodule generated by
   the labels.  In V (+) M' (M' with a large lowest-weight space near 3h) the modes preserve each summand:
     - labels in V: the ranks are those of V alone; M' is never reached, and the line theorem then forces dim V_{<=4h} >= N^2
       inside V itself, so V itself must have the ramp-plateau-ramp profile;
     - labels in the lowest-weight space of M': the budget restarts at the bottom of M' (E0 = 0, B = H) and (G) gives rho >= 3/4.
   The profile of a direct sum (a maximum of concave functions) is therefore not the profile the pieces see.  [P]
 * Exact ranks with a degenerate ground space follow the STEP profile, not a plateau [E, ramond_so8.out]: weights 0,1,2 have
   dimensions 16, 128, 576;  pieces (rr, rc) = (16, 448), (128, 128), (448, 16) for the 28 x 16 crossing (N N' = 448) and
   (16, 224), (128, 128), (224, 16) for a square N = 16 restriction (N^2 = 256): worst piece 4 + 7.8 = 11.8 of 16 bits,
   rho = 0.74, i.e. the 3/4 of (G) and of the reviewer's point (3) ("a pure step function costs 3/4").  Thresholds (16, 112), (112, 16).
 * What a plateau would need: ONE cyclic module with dim (V.W)_{<=E} = N^{1+o(1)} for E0 <= E <= E0 + H while N dual pairs of
   locality order H act on it.  The translation operator only forces dim M_{E0+k} >= N (polynomial growth, compatible).
   Not constructed, not excluded; it is equivalent in strength to P(N) = N^{1+o(1)}, so no obstruction theorem is expected
   from these methods.  [O]
 * Generators of weight >= 2 (reviewer's other non-concave case): the non-concavity is at weights below the first generator
   and goes the wrong way; by (L) a gain needs d/db[b log Z] > 0 somewhere, and every boson/fermion tower of weight
   1/2 <= D <= 8 has it strictly negative on the whole grid b in [0.01, 20] (locexp2.out) [N; a grid, not a proof].  Agrees with G20c (F3).
