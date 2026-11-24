# G23a — group-invariant (orbit) encoding of ST intercalate arrays  (working notes)
Tags: [E] exact computation, [P] proof (unreviewed), [H] heuristic, [O] open.
Model (g22b Thm 1): N x N array, cells coloured+signed, colours split S|T, every class a partial permutation,
every 2x2 with a mono diagonal has a mono other diagonal of the OTHER palette and sign product -1.  p = #T, q = #S.

## 1. Orbit encoding (orbit.py)
Gamma of order N acts regularly on rows and columns by left multiplication; cell (a, a c) gets colour (i(c), a t_c H_i):
block i(c) (palette fixed per block), H_i a NORMAL subgroup (colour stabiliser; H = 1 "free block", N colours per block,
H of order 2: N/2 colours), t_c in Gamma (canonical coset rep).  Class of colour (i, gH) = { (g h t_c^{-1}, g h t_c^{-1} c) }.
Q(c,d,v) := [cells (1,c) and (v, v d) have the same colour]  (N^3 vars, Q(c,d,v) = Q(d,c,v^-1)).
Theorem-1 conditions, reduced by invariance to x = 1:  Q(c,d,1) = 0 (c != d), Q(c,d,v) = 0 if v d = c (partial perm),
Q(c,d,v) -> Q(vd, v^-1 c, v)  and  pal(c) != pal(vd)   (other diagonal (1,vd),(v,c) mono, opposite palette).
Q defined both ways from t, i: ~ N^4 |H| + N^4 clauses (N=16: 2e5, N=32: 2e6) instead of ~N^4 * colours direct.
Signs: F_2 linear system on the full array after each colouring; blocking clause if inconsistent.
Every SAT output is re-verified: direct colour check, g22b shape.verify_shape (full [N]^4 identity + flint ranks),
g22c ev.py (p = rank(F - sum_S K (x) K^T) exactly, q = #terms; dense check at N <= 16).
Symmetry breaking: t_1 = 1 (valid for normal H), identical blocks ordered by min element; block caps |D|*|H| <= g*(other).
Validation [E]: order 8: only Q8 SAT for (8,8) (Z8, Z4xZ2, Z2^3, D4 UNSAT) = g22b; N=16 free 2+2 -> (32,32) SAT
(Q8xZ2, Pauli), exact via verify_shape and ev.py (dense) = 32/32.

## 2. The Q8 point in closed form  [E: q8closed.py, all 128 solutions enumerated by q8an.py]
Rows = cols = Q8 = {+-1, +-i, +-j, +-k}.  Cell (a,b), c = a^-1 b = eps*u, eps = +-1, u in {1,i,j,k}.
  sign(a,b) = eps                                    (left-invariant sign; found by solving the F_2 system restricted to s = sigma(c))
  palette   = S  iff c in {-1, i, j, k}  (a transversal of Q8/{+-1}),  T iff c in {1, -i, -j, -k}
  colour    = a * phi^-1(u),  phi = the order-3 automorphism i -> j -> k -> i  (index 0..7 within each palette)
Verified: colour conditions 0 violations, verify_shape identity True rr(R)=8 rc(C)=8, ev.py p=8 (dense 8), 16 classes all of size 4.
All solutions of the Q8 ansatz (t_1 = 1): 128 = 16 x 8; the S-block is ANY of the 16 transversals of Q8/Z(Q8), T-block = -S,
the remaining factor 8 is the colour relabelling of the second block.  In every solution t(c) = t(-c) (t factors through Q8/Z).
Structure: mod Z, Q8/Z = F_4^+ and phi = multiplication by omega; the colour index mod Z is a + omega(b - a) = omega^2 a + omega b,
a LATIN square over F_4; the central lift splits each mod-Z cell pair (c, -c) between S and T ("quaternionic doubling"
[[X, Y],[Y', X']]).  Zero slack: N^2 = 64 = q g*(p) + p g*(q) = 8*4 + 8*4; each T class (h = 4) spans a 4x4 submatrix with
6 = C(4,2) distinct S colours on symmetric pairs (of 8), and every pair of same-coloured cells spans an intercalate (48 intercalates).
Why it does not scale as-is: the mod-Z Latin-square structure forces class size |Q8/Z| = N/2, and N/2 <= g*(N) only for N <= 8.
The right generalisation for the 4/3 law is h = N^{2/3}: q = 2 -> 8 = 2^3, h = 4 = 2^2: Q8 = Sylow-2 of SU(3,2) (order q^3,
centre F_q, quotient F_{q^2}).  Next member: q = 4, Sylow-2 of SU(3,4), order 64, centre of order 4 -> 4 blocks = 4 transversals
of Z, 2 S + 2 T, each of 16 = q^2 elements, class size 16 = g*(128) -> would give exactly (128,128)@64 (the model lower bound).

## 3. Results, N = 16 (target below the direct sum (32,32); model bound 21)
Free colour orbits give p, q in 16 Z, so (24,24) needs stabilisers.  All UNSAT [E] (each < 75 s), all 14 groups of order 16
(abelian ones: Z16, Z4^2, Z8xZ2, Z2^2:Z4 complete; Z4xZ2^2, Z2^4 not run for f+h2):
  f+h2  (S: free + H2, T: free + H2, H2 central of order 2, all <= 12 unordered pairs (H_S,H_T)): UNSAT for every group run.
  h2x3  (S: 3 x H2, T: 3 x H2): UNSAT (Q8xZ2, Pauli, Q16, D4xZ2).
  f+h4h4 (S, T: free + 2 x H4, H4 normal of order 4): UNSAT (Q8xZ2, Pauli, Q16, D4xZ2, Z4:Z4, M16, Z2^2:Z4).
  semiregular Q8 x [2] (orbit2.py, Q8 acting with 2 row/col orbits, free colours, 3+3 blocks -> (24,24)): TIMEOUT 100 s;
  (16,16) UNSAT (consistent with the bound 21).
  (32,32) is SAT for Q8xZ2, Pauli, Z2^4, D4xZ2 (free 2+2), exact-verified (verify_shape + ev dense = 32/32).
  Free 2+2 (32,32)@16: SAT + exact-verified for 13 of 14 groups (incl. Z16, Z4^2, D8, Q16, ...); Z2^4: TIMEOUT 60 s.
  => least balanced value in the regular-group ansatz at N = 16: 32 for every group (= direct sum); nothing below found.

## 4. ZT ansatz ("Z-transversal", zt.py) = the Q8 structure generalised
Gamma, central Z of order m, Z0 < Z index 2. Blocks labelled by Z (palette S iff in Z0); every Z-coset meets every block once:
c = rep(cb) w lies in block lam(cb)^-1 w; t constant on Z-cosets.  Unknowns: tb (N/m cosets x N), lam (N/m x m) -> tiny SAT;
reduction (Gamma- and Z-invariance) makes every condition a condition on (coset pair, v).  m = 2: (N,N)@N; m = 4: (2N,2N)@N.
  N=8, m=2: Q8 SAT (reproduces the point, exact), D4/Z8/Z4xZ2/Z2^3 UNSAT.  N=16, m=4: (32,32) SAT for Q8xZ2 (exact).
  N=32, m=4 -> (64,64) (would beat the direct sum 128 and N^1.414 = 134): UNSAT for every (Z,Z0) tried (each < 7 s) in
  Q8xZ4, Q8xZ2^2, Q16xZ2, D4xZ4, Z4^2:Z2, Z8:Z4, Z4:Z8, Pauli x Z2, Z8xZ2:Z2, Z4^2xZ2, Z8xZ4, Z2^5  [E]
  (ES+/ES-, Q32, D16, SD32 have centre of order 2: ansatz not applicable with m = 4).
  N=64 = q^3, q = 4, Gamma = Sylow-2 of SU(3,4) (elements (a,b) in F16^2, b + b^4 = a^5, (a,b)(a',b') = (a+a', b+b'+a a'^4);
  q = 2 gives Q8 and the ansatz is SAT there for both mu): sub-ansatz t(cb) = (mu*a, *) (mod-Z part linear):
  UNSAT for all 14 mu in F16 \ {0,1} [E, 0.2 s each]; Z0 is WLOG (scalings act transitively on the 3 choices).
  Without the palette condition the same sub-ansatz IS satisfiable (but then same-palette intercalates appear, signs fail).

## 5. OBSTRUCTION for the linear ZT ansatz  [P, short proof, unreviewed; matches every SAT outcome]
Setting: 1 -> Z -> Gamma -> V -> 1 central, V = F_2^n, Z elementary abelian (additive), section s with s(0)=1 and bilinear
cocycle f (s(x)s(y) = s(x+y) f(x,y); such a representative exists for these extensions), so s(x)^2 = f(x,x) = "square map".
ZT data: c = s(x)(b + lam(x)) lies in block b in Z; palette S iff psi(b) = 0 (psi: Z -> F_2 linear, ker psi = Z0);
t(s(x)w) = s(tau x) sigma(x) with tau F_2-LINEAR on V (tau, tau+1 invertible are forced: rows/cols distinct), sigma, lam arbitrary.
CLAIM: a solution forces psi(f(a,a)) = 1 for every a != 0, i.e. the quadratic form a -> psi(s(a)^2) is ANISOTROPIC on V;
hence n <= 2 (Chevalley-Warning: every quadratic form in >= 3 variables over F_2 has a nontrivial zero).
PROOF. Take c, d in the same block b over cosets x = y + w, y (w != 0; ZT puts every coset in every block, so all pairs occur).
v = t(x)t(y)^-1 = s(tau w) A with A = f(tau w, tau y) + sigma(y+w) + sigma(y) (bilinearity); e = v d = s(tau w + y)(f(tau w,y) + A
+ b + lam(y)), whose block is f(tau w, y) + A + b + lam(y) + lam(tau w + y). Palette flip (Theorem 1) <=> psi(block(e) - b) = 1:
   E(w,y):  L_w(y) + g(y) + g(y+w) + Lam(y) + Lam(y + tau w) = 1   for all y and all w != 0,
 with L_w(y) = psi f(tau w, (tau+1) y) (bilinear in (w,y)), g = psi o sigma, Lam = psi o lam.  In +-1 notation (G = (-1)^g, ...):
   D_w G(y) D_{tau w} H(y) = -chi_w(y),  D_w G(y) := G(y)G(y+w),  chi_w(y) = (-1)^{L_w(y)} (a character in w and in y).
 For w1 != w2 both != 0: expand E(w1+w2,y) with D_{w1+w2}G(y) = D_{w1}G(y) D_{w2}G(y+w1) (same for H, tau linear), divide by
 E(w1,y) and by E(w2,y+w1):  H(z)H(z+a)H(z+delta)H(z+a+delta) = -(-1)^{psi f(a,delta)} for all z, with a = tau w2,
 delta = (tau+1) w1.  Choose w2 = tau^-1 a, w1 = (tau+1)^-1 a (nonzero, distinct for a != 0): then delta = a and the left side
 is H(z)^2 H(z+a)^2 = 1, so psi f(a,a) = 1 for all a != 0.  QED
 (Even without the a = delta trick: Lam has constant second derivatives, so Lam is quadratic with polar form 1 + psi f on almost
 all pairs, which is impossible for n >= 3 by counting values of a linear functional.)
Interpretation: every non-central element must square into the T-part Z \ Z0.  Q8: all 6 non-central squares = -1, Z0 = {1}. OK.
Checks against SAT [E]: Q8xZ2 (Z = <-1> x Z2, n = 2): squares of non-central elements = (-1,0) =: element 8; the ZT runs gave
SAT for Z0 = {0,1} and {0,9} and UNSAT exactly for Z0 = {0,8}. N = 32, |Z| = 4 (n = 3): all UNSAT (also for nonlinear tau).
SU(3,4) (n = 4): squares (0, a^5) hit all of Z \ 0 -> isotropic -> UNSAT for all 14 mu, as observed.
CONSEQUENCE: the ZT family (the literal generalisation of the Q8 point) has n <= 2, i.e. blocks of size |V| <= 4: classes of
size <= 4 and p ~ N^2/8 -- no better than tiling.  The linear-tau / bilinear-cocycle hypotheses are where the proof needs
them; the general (nonlinear tau) ZT SAT at N = 32 is UNSAT for all groups tried as well [E].
[E, anis.py] all (Z, Z0) of all groups of order 8 (m=2) and 16 (m=4), GENERAL tau: every ZT-SAT case satisfies the
anisotropy criterion "non-central squares avoid Z0" (Q8; Q8xZ2 with Z0 = {0,1},{0,9}); it is not sufficient
(Z8, Z16, Z8xZ2, Z4:Z4 anisotropic but UNSAT).  So the criterion is empirically necessary beyond the linear-tau proof.
## 3b. N = 16, all counting-feasible stabiliser specs with max(p,q) <= 24  (all16.py; free / H2 / H4 blocks, first normal
H2 and H4 of each group, all 22 palette specs passing the Key-Lemma capacity count) [E]:
  Q8xZ2, Pauli, Q16, D4xZ2: every spec UNSAT (84 of 88 runs, < 20 s each), except f+h2 for Q8xZ2/Q16 = TIMEOUT 30 s under load,
  but f+h2 was proved UNSAT for these groups (all H-pairs) in sect. 3.  => no regular-group array with max(p,q) <= 24 at N = 16
  in this (normal-stabiliser) ansatz for these 4 groups; other groups covered only by the f+h2 / h2x3 / f+h4h4 families.

## 6. Table (least balanced max(p,q) found in the orbit ansatz; model lower bound from Thm 3; direct sum N^2/8)
 N  | model LB | direct sum | N^1.414 | result
 8  |   8      |    8       |   19    | (8,8) Q8 only (ZT and free); Z8, Z4xZ2, Z2^3, D4 UNSAT
 16 |  21      |   32       |   50    | 32 for all 14 groups (free 2+2; Z2^4 timeout); <= 24: UNSAT in every spec run (sects 3, 3b)
 24 |  36      |   72       |   90    | semiregular Q8 x [3], (40,40): see sect. 7
 27 |  41      |  (~91)     |  106    | free 2+2 (54,54): Z9:Z3 UNSAT (14 s); Heis3, Z3^3, Z9xZ3: TIMEOUT 150-500 s (undecided)
 32 |  52      |  128       |  134    | (64,64): ZT ansatz UNSAT for all 12 groups with |Z| >= 4; general free 2+2 (2.1M clauses):
    |          |            |         | ES-, ES+, Q16xZ2 TIMEOUT (150-580 s, with and without the z-doubling symmetry)
 64 | 128      |  512       |  358    | (128,128): SU(3,4)-Sylow linear-tau ZT UNSAT (all mu), proved impossible (sect. 5);
    |          |            |         | general ZT (5M clauses): see sect. 7; Q8xQ8 ZT TIMEOUT 140 s
NB at N <= 34, N^2/8 < N^1.414 numerically: beating "N^1.414" at a single small N says nothing about the exponent.

## 7. More N = 27 / 32 / 64 runs  [E]
 * (64,64)@32 general free 2+2 (2.1M clauses): ES- TIMEOUT 280 s; ES+ TIMEOUT; with z-doubling (t(zc)=t(c), block(zc)=partner):
   ES-, Q16xZ2 TIMEOUT 150 s; ES-, ES+ with z-doubling AND exact block sizes 8: TIMEOUT (~9 min, killed).  Undecided.
 * Same with equivariance under an automorphism phi (t(phi c) = phi t(c), blocks phi-invariant; the Q8 point IS equivariant
   under its order-3 automorphism): ES- order 5 (3 automorphisms, with/without z-doubling) UNSAT 2-4 s; ES+ order 3 (z-doubling)
   UNSAT.  Order-2/3/4 runs on many groups: see A32o*.txt (summary below).
 * (54,54)@27 Heis3: plain TIMEOUT 500 s; with automorphism equivariance of order 2, 3, 4, 8 (3 automorphisms each): all UNSAT (< 8 s).
   Z9:Z3 plain UNSAT 14 s; Z3^3 TIMEOUT 150 s.
 * (128,128)@64 SU(3,4)-Sylow, general (nonlinear tau) ZT (5.0M clauses): TIMEOUT 560 s.  Q8xQ8 ZT: TIMEOUT 140 s.

## 8. Verdict
 * Logic check: any ST array (Thm 1 of g22b) IS an exact identity W^T = sum_S K W K + sum_T tr(KW) K with q = #S terms and
   rank R = #T = p; ev.py recomputes p = rank(F - sum_S K (x) K^T) over Q independently of the model.  So a family of arrays with
   max(p,q) = O(N^beta), beta < sqrt2, WOULD be a genuine record for the general problem.  None was found.
 * Q8 point: closed form above (sect. 2), = Sylow-2 of SU(3,2) with the omega-scaling; zero slack.
 * Group (orbit) ansatz: no point below the direct sum anywhere: N=16 best 32 (all groups), every (24,24)-type stabiliser spec
   UNSAT; N=32 (64,64): UNSAT in the ZT ansatz (all groups with |Z| >= 4) and in every automorphism-symmetric variant tried,
   plain instances undecided (timeouts); N=27 (54,54): UNSAT (Z9:Z3; Heis3 with any automorphism symmetry), plain Heis3 undecided;
   N=64 (128,128) in SU(3,4): linear-tau ZT impossible (proof), general ZT undecided.
 * Obstruction [P]: the literal generalisation of the Q8 point (ZT, linear tau) requires psi(x^2) != 0 on Gamma/Z \ 0, an
   anisotropic quadratic form, hence rank(Gamma/Z) <= 2: block size <= 4, no family.  The SU(3,q) / q^3 family suggested by the
   4/3 counting is killed at q = 4 for exactly this reason (Chevalley-Warning).
 * Open: plain free 2+2 at N = 32 (ES+-, Q16xZ2, ...), N = 27 Heis3, general-tau ZT at N = 64; generic (non-group) ST arrays
   below N^2/8 at any N >= 16 (g22b: (28,28)@16 undecided).

## 9. Final run log (supersedes "see A32o*" in sect. 7)  [E]
 * Automorphism-equivariant free 2+2 (64,64)@32, 3 random automorphisms per order, 50 s limit: ALL UNSAT
   order 2: D16, ES+, ES-, Q16xZ2, Q32, Q8xZ2^2, Q8xZ4, SD32;  order 4: same + Q8oZ8 (=Z8xZ2:Z2), Z4:Z8, Z4^2:Z2, Z8:Z4;
   order 3: Pauli x Z2, Q8xZ2^2, Q8xZ4, Z2^5, Z4^2xZ2;  order 5: ES- (also with z-doubling).
 * Plain ES- free 2+2 with kissat404: TIMEOUT 580 s.  Plain Heis3 (54,54)@27 with kissat404: TIMEOUT 580 s.
 * Semiregular Q8 x [3] (40,40)@24: TIMEOUT 580 s.
All processes killed at the end.
