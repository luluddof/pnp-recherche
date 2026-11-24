# G24g — strong swap lemma, small-q regime (working notes, saved as I go)
Tags: [P] proof here (short, unreviewed), [E] exact computation, [N] numerical, [H] heuristic, [O] open.
f(N,r) := min rank(Swap - K) over K of operator-Schmidt rank <= r  (rank in (x,y)|(u,v)).

## 1. Reformulation used throughout
Swap - K = (1 - Theta) Swap, Theta = K Swap, so rank(Swap - K) = rank(1 - Theta).  In coordinates
(V_{xy} <-> e_x (x) e_y): (Theta V) has row b equal to sum_{b'} sum_l (E_l)_{b b'} C_l (column b' of V), i.e.
  Theta[(b,j),(k,b')] = sum_l (E_l)_{b b'} (C_l)_{j k}.
Symmetry group [P]: for X, Y invertible, M -> (X(x)Y) M (Y^{-1}(x)X^{-1}) fixes Swap and sends
E(x)C -> X E Y^{-1} (x) Y C X^{-1}; rank and osr are preserved.  So only the E-pencil up to strict equivalence matters
for the argument below.
Principle: rank(1 - Theta) >= rank of any principal compression (1 - Theta)[Om,Om] = |Om| if Theta[Om,Om] is nilpotent.

## 2. TRIANGULAR SWAP LEMMA  [P, any field, any characteristic; exact checks tri_chk.py]
Lemma T.  If there are invertible X, Y with X E_l Y^{-1} upper triangular for all l, then for ALL C_1..C_r
            rank(Swap - sum_l E_l (x) C_l) >= N(N-1)/2.
Proof. After the symmetry, E_l upper triangular.  Om = {(b,j) : b > j} (|Om| = N(N-1)/2).  For (b,j), (k,b') in Om the
entry Theta[(b,j),(k,b')] = sum_l (E_l)_{b b'} (C_l)_{jk} is nonzero only if b' >= b > j.  Order Om by the second
coordinate, decreasing: Theta[Om,Om] is strictly triangular (an entry goes from column-index b' to column-index j < b'),
hence nilpotent, (1-Theta)[Om,Om] is unitriangular, rank >= N(N-1)/2.  QED
Lemma P (pencils).  Over an algebraically closed field every pair (E1,E2) of N x N matrices is strictly equivalent to
a pair of upper triangular matrices.  Proof: Kronecker canonical form = block diagonal with blocks L_eps (eps x (eps+1),
[I|0],[0|I]), L_eta^T, Jordan blocks (I, J_k(lam)) and (N_k, I).  Order the blocks: all L_eps, then the regular blocks,
then all L_eta^T (as many L^T as L since the pencil is square).  With row offset p and column offset q of a block
(natural order inside), L_eps needs p <= q and raises q - p by 1, regular blocks need p <= q, L_eta^T needs p + 1 <= q
and lowers q - p by 1.  Starting from p = q = 0 all conditions hold.  QED  (checked on 375 random block lists, tri_chk.py)
THEOREM 1 [P].  For every field K and every N:  f(N,2) = f(N,1) = N(N-1)/2.  i.e.
            rank(Swap - E1(x)C1 - E2(x)C2) >= N(N-1)/2,  tight at E1 = C1 = I, E2 = C2 = 0.
(pass to the algebraic closure: ranks do not change.)  r = 1 is the special case of one matrix (rank normal form).
 * This removes the char != 2 hypothesis of G19b Thm 3 (r = 1), and settles r = 2 (G19b had only rank >= N, Thm 4).
 * Consistent with: G24d over F_2: (5,2)@4 UNSAT, (6,2)@4 SAT; numerics N = 4, 5 (r = 2 -> N(N-1)/2).
 * Not a limit statement: border (approximate) decompositions are not covered (as for G19b Thms 3, 4).
   Border (p,2) means Swap = lim (A_t + B_t) with rank A_t <= p, osr B_t <= 2, where A_t, B_t may diverge separately;
   the proof says nothing about such limits.  So border (p,2) with p < N(N-1)/2 is NOT excluded.  [O for border;
   numerics (G19b): no border gain at q = 2 for N = 4, 5]
Staircase generalisation (Lemma C) [P]: for arbitrary E_l, choose bases; let t(b) = number of leading columns in which row b
of every X E_l Y^{-1} vanishes.  Om' = {(b,j): j < b, j < t(b)} is nilpotent by the same argument, so
            rank(Swap - K) >= max_{X,Y} sum_b min(b-1, t(b))      (b = 1..N).
Since dim (sum_l E_l F) <= r dim F, one can always reach t(b) >= floor((b-1)/r) (flag F_t, rows in (E F_t)^perp):
            f(N,r) >= sum_{b=1}^N floor((b-1)/r)  ~  N^2/(2r) - N/2       [P]
(Theorem C gave max_k k(N-kr) ~ N^2/(4r): a factor 2 better, any field.)  CAVEAT: this is E-side only; for an
E-span that is a dimension expander no flag does better than ~N^2/(2(1+eps)), so Lemma C cannot reach N^2/2 - crN
for r >= 3 by itself.
 IMPROVEMENT with the pencil [P]: by Lemma P two elements E1, E2 of the span can be made upper triangular, so for the flag
 F_t = first t columns, dim(E1 F_t + E2 F_t) <= t and dim(sum_l E_l F_t) <= (r-1) t.  Hence for r >= 2
            f(N,r) >= sum_{b=1}^N floor((b-1)/(r-1))   (r = 2: N(N-1)/2 = Theorem 1; r = 3: ~ N^2/4).   any field.
 Exponent-wise this is the same as Theorem C (N^2/r up to constants).

## 3. Where exactly r >= 3 escapes  [P]
 * Lemma T and its mirror (conjugate by Swap: E(x)C -> C(x)E) show: an exact gain below N(N-1)/2 needs BOTH spans
   span{E_l} and span{C_l} to be non-triangularizable by equivalence.  Dimension count: an r-dim space of N x N matrices
   is generically triangularizable iff r <= 2.  So r = 3 is the first case where a gain is possible, consistent with
   f(4,3) = 5 (G19b exact point).
 * Canonical form for a generic 3-dim span [P]: (E1,E2) ~ (I,D); the pairs (X,Y) = (b P_s, b' P_s) (b, b' in the Borel B,
   s a permutation) keep them upper triangular, and E3 -> b M b'^{-1}: by the Bruhat decomposition E3 becomes a permutation
   matrix w (generically the antidiagonal w0).  So generic span{E} ~ span(U1, U2, w0); Lemma C then gives t(b) = min(b-1,N-b),
   sum ~ N^2/4: the E-side alone cannot do better (span(I,D,w0) is far from triangularizable), the C-side must be used.
 * The exact (5,3)@4 point does NOT contain a term G(x)G^{-1} with coefficient 1 (emb53.py: min residual of G^{-1} off
   span(G') = 0.5), so it does not embed into all N as "I(x)I + osr-2 correction"; I found no construction giving an exact
   gain at r = 3 for N >= 5.

## 4. Data  (least rank f(N,r); exact = verified decomposition + proof, border = numerical limit only)
 r = 1: N(N-1)/2 [P, all fields now].   r = 2: N(N-1)/2 for ALL N [P, Theorem 1] -- "no gain at all over r = 1",
 exactly as the N = 4, 5 numerics suggested (task question answered: yes).
 r = 3: N = 3: 1 (mirror of (3,1));  N = 4: 5 exact (G19b);  N = 5: exact <= 10, BORDER 8 (G19b: (8,3), (9,3) tails 3e-7
 with coefficients 60-80; (7,3) tail 2e-3 unresolved); lower bound proved: 4 (Lemma C+pencil, sum floor((b-1)/2)).
 N = 6, 7: runs below (vp.py copy g24gswp.py, machine load ~20 on 4 cores, so few restarts).
 Fit question: with r = 2 flat, the exact data are compatible with f = N(N-1)/2 - O(r) (gains 0,1,2,3 at N = 4,
 r = 2..5; N = 5: r = 4 -> 7 exact); the border data at N = 5 (r = 3 -> 8) show a gain of 2 already at r = 3.
 No fit of the form N(N-1)/2 - alpha(r) N + beta is supported: alpha(2) = 0 exactly, and for r >= 3 there are too few exact
 points (N = 4 only) to separate alpha(r) N from a constant.

## 5. Verdict and consequences
 PROVED [P, any field]:  (i) Lemma T: E-family (or C-family) simultaneously triangularizable by equivalence
   => rank(Swap - sum_l E_l(x)C_l) >= N(N-1)/2, for EVERY r.  (ii) Theorem 1: f(N,2) = f(N,1) = N(N-1)/2 (r = 1 now also
   in characteristic 2).  (iii) f(N,r) >= sum_{b<=N} floor((b-1)/(r-1)) ~ N^2/(2(r-1)) (factor ~2 over Theorem C).
 NOT PROVED: the strong lemma f(N,r) >= N^2/2 - c r N for r >= 3.  Obstruction located precisely: for r >= 3 both spans
   can be non-triangularizable (generic span{E} ~ span(U1, U2, w0) by Bruhat), the E-side staircase then stops at ~N^2/4,
   and any proof must use the coupling between the E- and C-sides.  No counterexample either: the exact r = 3 gains
   found (sect. 7-8: 3 at N = 6 and N = 9, from A_3 (x) I_m) are bounded; linear gains exist at r = 4 (N/2, Swap_2 (x) I),
   which the strong lemma allows.  [Section 5 written before sect. 7-10; see those for the final data.]
 G24f window (E0 = 2, 3): NOT closed by this.  The general bound (iii) has the same exponent as Theorem C
   (N^2 <~ r D, kappa <= E0 - 1/2), and dim M_{<E0} >= 1 + 3n > 2 there, so Theorem 1 does not apply as such.
   BUT Lemma T gives a checkable criterion with no restriction on r: if G24f's actual family {E_m : m in M_{<E0}}
   (P_b = sum_m E_m(x,y) C_m(v,u)) - or the C_m family - is simultaneously triangularizable by equivalence (e.g. because
   the E_m respect a grading/filtration of the index spaces with all shifts on one side), then N(N-1)/2 <= dim M_{<=E0}
   for every E0, which is exactly the E0 = 1 conclusion (kappa <= E0/2) and would close the whole H = 0 window.
   Whether their family has this property I did not check (needs their build.py objects).  [O]
 Trade-off curve near the endpoint (exact frontier, symmetric under p <-> q): flat for two steps at both ends:
   (N(N-1)/2, 1), (N(N-1)/2, 2) and mirror (1, N(N-1)/2), (2, N(N-1)/2) [P]; the first possible drop is at q = 3
   (realised at N = 4: (5,3)).  Border frontier: Theorem 1 does not cover limits; numerics show no border gain at q = 2
   ((5,2)@4 tail 1.05, (9,2)@5 tail 1.07, G19b) but a border gain at q = 3 ((8,3)@5).
- GF(2), N = 3, r = 2 [E, gf2_r2_n3.py]: exhaustive over 35 pencil orbit reps (E1 = diag(I_k,0)) x all 2^18 (C1,C2): min rank = 3 = N(N-1)/2 (independent check of Theorem 1 in characteristic 2).

## 6. Lemma C in flag form (the usable certificate)  [P, any field]
 For any complete flag 0 = F_0 < F_1 < ... < F_N = K^N (dim F_t = t) and E-span 𝓔 = span{E_l}:
        rank(Swap - sum_l E_l (x) C_l) >= sum_{t=1}^{N-1} ( N - max(t, dim 𝓔F_t) )
                                         = N(N-1)/2 - sum_{t=1}^{N-1} (dim 𝓔F_t - t)_+ ,
 for ALL C_l; the same with 𝓒 = span{C_l} (conjugate by Swap: Swap K Swap = sum C_l (x) E_l).
 Proof: rows x_b chosen in an adapted basis of the decreasing chain (𝓔F_t)^perp, bottom rows deepest; then
 #{b : min(b-1, t(b)) >= t} = min(N - dim 𝓔F_t, N - t) and Lemma C's nilpotent set has that many elements per level.
 So the gain below N(N-1)/2 is at most the "flag defect" min_F sum_t (dim 𝓔F_t - t)_+ (and the same for 𝓒).
 Strong lemma <= any family with flag defect <= c r N on one side.  False as a route for dimension-expanding spans
 (defect ~ N^2/4 for generic r = 3), but it is a COMPUTABLE certificate for a concrete family:
 e.g. if E_m = A O_m B with O_m all preserving one flag of a module M (A: M -> K^N, B: K^N -> M), take F_t from the
 preimages B^{-1}(M_i) and bound dim 𝓔F_t by dim A(M_i).  Suggested check for G24f's actual P_b family.

## 7. NEW EXACT POINT: r = 3 gains 3 at N = 6  [E, verify_exact.py, exact_6_12_3.json]
 f(6,3) <= 12 = N(N-1)/2 - 3.  Found by vp (complex: 1/3 restarts converge in 83-275 its with coefficients 1.3; real
 version vpreal.py: 3/8 restarts), gauge-sparsified (sparsify.py, copy of G19b's) to a MONOMIAL +-1 solution, rounded,
 verified with flint: rank_Q(Swap - K) = 12, osr_Q(K) = 3; same ranks mod 2, 3, 5, 7 (so valid over every field).
 Convention K[x,y,u,v] = sum_l E_l[x,u] C_l[y,v], Swap[x,y,u,v] = [x=v][y=u]; E_ij = matrix unit (1-based):
   E_1 = -E23 +E35 +E41 -E52    C_1 = +E14 -E25 -E32 +E53
   E_2 = -E15 -E44 +E56 -E63    C_2 = -E36 -E44 -E51 +E65
   E_3 = -E11 +E26 +E34 +E62    C_3 = -E11 +E26 +E43 +E62
 (each E_l, C_l a signed partial permutation of rank 4.)  Spectrum of Theta = K Swap (numerical solution):
 1^24, (-1)^8, (-2)^3, 2^1.
 => r = 3 is NOT "N(N-1)/2 - O(1) only at N = 4": the gain at r = 3 is 1 (N=4), >= 0 exact / 2 border (N=5), >= 3 (N=6).
 Remarks: (a) the ST (monomial, no-cancellation) model of G24d says (12,3)@6 UNSAT (st/o_6_12_3.txt, all s): the new point
 is monomial but uses CANCELLATIONS between R = Swap - K and K (supp K not inside supp Swap), so it lies outside the ST
 model; ST lower bounds are not lower bounds for f.  ST also gives (14,2)@6 UNSAT, consistent with Theorem 1.
 (b) vp scan (complex, OMP_NUM_THREADS=1, 3000 its; logs/vp_*_s11.txt): N=6: p=11 d2 6e-5 (coef 4-10, still decreasing,
 unresolved), p=10 d2 0.19, p=9 0.40 (not found); N=5: p=9,8,7 d2 3e-4, 1.4e-3, 5e-3 (coef 10-40: border-type, as G19b);
 N=7: p=18 d2 4e-3, p=16 1.3e-2 (unresolved, 4 restarts).

## 8. The (12,3)@6 point explained; tensor families  [E, tensor_fam.py; P for the rank formula]
 Phi = sum_x e_x (x) e_x.  A_k := Swap_k - Phi Phi^T has osr k(k-1)/2 (realignment: Swap - I(x)I) and rank(Swap_k - A_k) = 1
 (this is the mirror (1, k(k-1)/2) of the r = 1 point).  On K^N = K^k (x) K^m, Swap_N = Swap_k (x) Swap_m and for
 K = A (x) (I_m (x) I_m) (osr = osr A):   rank(Swap_N - K) = m(m+1)/2 * rank(Swap_k - A) + m(m-1)/2 * rank(Swap_k + A).
 * A = A_3, m = 2: 3*1 + 1*9 = 12 at N = 6, osr 3  -- the same value as the vp/sparsified point (which has the same
   Theta-spectrum structure); exact check tensor_fam.py: rank 12, osr 3.  m = 3: 33 at N = 9 (gain 3); m = 4: 64 (gain 2);
   m = 5: no gain.  So this r = 3 family has BOUNDED gain (3).
 * A = Swap_k (osr k^2): rank = k^2 m(m-1)/2 = N(N-k)/2, i.e.  f(N, k^2) <= N(N-1)/2 - (k-1)N/2  for k | N  [E: k=2,3,4].
   In particular r = 4: gain N/2 (N even) -- a LINEAR gain at constant r: the strong lemma cannot have c < 1/8 at r = 4,
   and in general gains of order sqrt(r) N/2 occur.  Consistent with N^2/2 - c r N (c = 1/2 suffices for this family).
 * Large-m criterion [P]: the family gives a gain linear in N iff x + y = k^2 with x := rank(Swap_k - A) < k(k-1)/2,
   y := rank(Swap_k + A); x + y = k^2 holds when Theta_A = Swap A is an involution, and then gain = m (k(k-1)/2 - x).
   So a linear-gain family at osr r <=> an osr-r A with (Swap A)^2 = 1 and tr(Swap A) = sum_l tr(E_l C_l) > k.
   (For r <= 2 impossible by Theorem 1; r = 4: A = Swap_2, tr = 4 > 2.)  Open for r = 3.

## 9. Numerics for N >= 7 (r = 3) and r = 4  [N; real solver vpreal.py, warm starts warm.py, continuation cont.py]
 Random starts never converge to exact points for N >= 7 at r = 3.  Best values (d2, max coefficient, all still slowly
 decreasing with slowly growing coefficients = border-type, unresolved):
   (11,3)@6 1.0e-5 (coef 12);  (18,3)@7 4.8e-4 (33);  (17,3)@7 1.2e-3 (32);  (25,3)@8 5.8e-5 (16);  (24,3)@8 1.3e-2;
   (27,3)@8 1.1e-3;  (22,3)@8 3.3e-2 (not found).  r = 4 at N = 6: (11,4) 2.1e-4, (10,4) 1.9e-2.
 Fixed E-span tests (fixE.py): with span{E} = span(I, D, w0) or random 3-dim, optimising only the C's gives NO gain at
 N = 4, 5 (tails 1.46 / 0.99 at p = 5, N = 4; 1.48 / 1.37 at p = 9, N = 5): gains need special (non-generic) E-spans.

## 10. Data table (least rank f(N,r); E = exact decomposition known, P = proved lower bound, b = border/numerical only)
   N | r=1  | r=2   | r=3                                   | r=4
   3 | 3 P  | 3 P   | 1 E                                   | 1
   4 | 6 P  | 6 P   | 5 E (P >= 2)                          | 4 E
   5 | 10 P | 10 P  | <= 10 E, b 8 (P >= 4)                 | 7 E
   6 | 15 P | 15 P  | <= 12 E (new), b 11? (P >= 6)         | <= 12 E (Swap_2 (x) I), b 11?
   7 | 21 P | 21 P  | <= 21 E, b 18? (P >= 9)               | -
   8 | 28 P | 28 P  | <= 28 E, b 25? (P >= 12)              | <= 24 E
   9 | 36 P | 36 P  | <= 33 E (A_3 (x) I_3)                 | -
 Fit: r = 2 is EXACTLY flat (alpha(2) = 0, beta = 0).  r = 3: exact gains 1, 0?, 3, ?, ?, 3 (N = 4..9), border gains
 2, 4?, 3?, 3?  -- bounded within N <= 9, no sign of alpha(3) N.  r = k^2: alpha >= (k-1)/2 (exact family).
 GF(2) exact (9,3)@5 (G24d gf2.py copy in gf/, GNF = 3 and 5, CryptoMiniSat): no answer within the 560 s limit.  [O]

## Files
 notes.md; tri_chk.py (Lemma T / Lemma P / staircase checks, exact); flag_chk.py (flag form); gf2_r2_n3.py (exhaustive
 GF(2), N = 3, r = 2); exact_6_12_3.json + verify_exact.py + sparsify.py + vpreal.py (new exact point); tensor_fam.py
 (A_k (x) I, Swap_k (x) I families); emb53.py; fixE.py; warm.py, cont.py, g24gswp.py (= vp.py copy); logs/; st/ (ST runs);
 gf/.  All background jobs killed at the end.
