# G18h — is the floor on rho an artefact of the "single graded space" framework?  (working notes, unreviewed)

Status tags: [P] proved here (short argument written), [N] exact computation, [Np] computation mod a prime
(not exact over Q), [H] heuristic / accounting only, [O] open.
Sources: Claude's paper (Def. 7, Thm A, App. B, App. H), G9a, G9b, G6b, N1 notes, G18c's "D^{-1/2} law" (relayed
by the coordinator, heuristic).  Nothing here is a new algorithm or a lower bound for 3-SAT.

BOTTOM LINE.  The floor rho >= 1/sqrt2 of G9a is a statement about an ACCOUNTING (ranks bounded by dimensions of
weight windows) plus one structural hypothesis taken from vertex algebras, (V).  As a restriction on pieces, "factor
through the graded components of one space" is EMPTY (Prop. 1): all the content is in (V) and in the shape of the
growth profile.  Dropping (V) or concavity, the rank count only gives rho >= 1/2, and explicit PROFILES reach 1/2 -
but I found NO operators realising them, and no piece family with asymptotic rho < 1/sqrt2.  The only floor proved
for all pieces stays rho >= 1/2.  The exponent of the algorithm is unchanged.

## 1. Definitions (restated)

Crossing: X(x,y,u,v) = [x=v][y=u], x,y,u,v in [N], N = 2^lambda.  Two flattenings: row (x,y)|(u,v), column
(x,u)|(y,v); X is a permutation matrix (rank N^2) in both.
Pieces (Def. 7, over a field F): for every lambda, functions A_i, i in I, with
 (a) sum_i A_i = d X, d != 0;  (b) log rank_row(A_i) + log rank_col(A_i) <= 4 rho lambda + O(log lambda);
 (c) |I| = poly(lambda), entries/factorisations computable in 2^{O(lambda)}.
rho = lim sup of max_i (c_row + c_col)/(4 lambda).  Theorem A: e_{l+1} = (rho + o(1)) e_l while K is squared, so
time 2^{O(m/(log m)^{log2(1/rho)})}: rho = 11/12 -> exponent 0.125; 1/sqrt2 -> 1/2; 1/2 -> 1 (m/log m).
General floor [P, N1/G9b]: rank is subadditive, so some piece has row rank >= N^2/|I|: rho >= 1/2 over every field.

The G9a framework ("graded window"), one assumption per line:
 A1 (one space)     a graded space W = (+) W_E with bases; labels y,u = N basis vectors of W_{E0}.
 A2 (operators)     x, v name operators a_x^(i)[r], b_v[r] : W_E -> W_{E+r}.
 A3 (both orders)   piece s is <u|Q_s(x,v)|y>, Q_s = sum_{i<=H} (-1)^i a^(i)[s-h+i] b[h-s-i] = same with the factors
                    exchanged (Lemma 9a);  sum_s Q_s = d [x=v] Id on W_{E0} (Lemma 9bc + 10).
 A4 (accounting)    the cost of a piece is bounded ONLY by dimensions: row-wise it passes through weights <= x_s,
                    column-wise through weights <= M - x_s, M = 2E0 + H ("what one order gains the other loses");
                    L(E) = log dim W_{<=E}.  (G18c [N]: for App. B the true ranks equal dim W_{<=x} until they
                    saturate at N^2, so the accounting is tight there.)
 A5 (interval)      the shifts s form an interval, every intermediate x_s is a piece that is charged.
 A6 (V)             state-field: N independent states of weight <= h_a give the operators, and H >= wt a + wt b;
                    with A1: E0 >= h1 and H >= 2 h1, where L(h1) = lambda.  Hence M >= 4 h1.
 A7 (shape)         L concave / power law c E^alpha / log L concave in log E.
 A8 (field)         Q (App. B needs d = h! invertible: characteristic 0 or > h).
 A9 (counting)      exact identity (a) over the ring where one counts; decision = a fixed finite field.
 A10 (poly pieces)  |I| = poly(lambda).
 A11 (uniform)      the same pieces at every level; D = 2 directions (rows, columns).

Floor proof of G9a, step by step, with what each step uses:
 S1 rank count: some piece has L(x_s) >= 2 lambda - log(|I|(H+1)), another has L(M - x_s') >= the same.
    Uses A4, A10 and subadditivity of rank (any field: A8, A9 not used).
 S2 balanced piece: by A5 and monotone L some charged piece has x ~ M/2 (or a piece is charged 4 lambda):
    cost >= min(4 lambda, 2 L(M/2)).  Uses A4, A5.
 S3 size of the budget: L(M) >= 2 lambda (S1) and M >= 4 h1 (A6).
 S4 shape: for L = c E^alpha, S1+S2 give rho >= 2^{-alpha}, S2+A6 give rho >= 2^{alpha-1}; max >= 1/sqrt2, equality
    only at alpha = 1/2.  Uses A7.  With merely concave L: 2/3 (G9a).
 A1-A3 are used only to JUSTIFY A4-A6 (they are true for App. B); they are not used as restrictions.

## 2. What each assumption is worth

### Prop. 1 [P]  A1-A3 alone are no restriction.
Let (A_s), s = 1..S, be ANY pieces with row factorisations A_s = sum_{tau<=r_s} f_tau(x,y) g_tau(u,v) and column
factorisations A_s = sum_{sigma<=c_s} f'_sigma(x,u) g'_sigma(y,v).  Put W_{E0} = F^N, W_{E0+s} = F^{r_s},
W_{E0-s} = F^{c_s}, and define (all other components zero)
  a_x[s]: W_{E0} -> W_{E0+s}, e_y -> sum_tau f_tau(x,y) e_tau;      b_v[-s]: W_{E0+s} -> W_{E0}, e_tau -> sum_u g_tau(u,v) e_u;
  b_v[-s]: W_{E0} -> W_{E0-s}, e_y -> sum_sigma g'_sigma(y,v) e_sigma;  a_x[s]: W_{E0-s} -> W_{E0}, e_sigma -> sum_u f'_sigma(x,u) e_u.
Then <u|b_v[-s] a_x[s]|y> = A_s = <u|a_x[s] b_v[-s]|y> (both orders agree on W_{E0}, with H = 0) and
sum_s a_x[s] b_v[-s] = d [x=v] Id.  So "pieces = matrix entries of a both-orders operator identity on one graded
space" (with H = 0 and an arbitrary dimension profile) is the general notion of pieces.  The floor 1/sqrt2 can
therefore only come from A6 (V) and A7 (shape), and it is a floor for vertex-algebra-like realisations only.

### Table of accounting floors [P for the inequalities; profiles checked numerically in section 4]
 kept assumptions                          floor of the accounting       extremal profile
 A4,A5,A6, power law                       1/sqrt2                       c sqrt(E) (partitions)
 A4,A5,A6, L concave                       2/3 (G9a)                     slope drops by 3 at h1
 A4,A5,A6, L only nondecreasing            1/2                           staircase: 0 | lambda on [h1,3h1] | 2 lambda on (3h1,4h1]
 A4,A5, power law, no (V)                  2^{-alpha} -> 1/2 at alpha=1  L linear (dim W_E = 2^E), E0 = lambda, H = o(E0)
 A4 only (no interval, no conservation)    1/2                           two pieces (N^2,1), (1,N^2) - the general floor
So in the accounting, (V) and the shape are each "proof-essential": remove either and the proof stops at 1/2.
The price: with (V), L must be flat (= lambda) on [h1, 3h1] and then double; without (V), an exponentially growing
graded space must carry N x N "almost commuting" (H = o(E0)) operator families.  I know neither.

### (a) other spaces / filtrations / non-abelian structures
- Filtration instead of grading: INESSENTIAL [P].  S1-S4 use only the nondecreasing function L(E) = log dim W_{<=E}
  and the conservation x_row + x_col = M.  App. B itself only uses the filtration (window = all weights <= x).
- Two spaces (row order through W, column order through W' with profiles L, L'): S1, S2 survive with
  L(x) + L'(M - x); the floor of the accounting is again 1/2 iff both profiles are linear.  The single space
  enters only through A6 (the labels and the operators live in the same space, so M >= 4 h1).  [P, accounting]
- Tensor products of piece families (two colours, rank-n Heisenberg): costs add, rho is the weighted mean; n-colour
  Fock spaces have Cardy growth exp(c sqrt(nE)) for E >> n: alpha = 1/2 again.  For E << n the growth is faster
  (alpha ~ 1), which is WORSE under (V) (2^{alpha-1} -> 1).  [P for the product; H for the growth statement]
- Clifford / octonions (App. H (74)): each product stands in ONE order, so a term has full rank in the other
  direction (rho >= 3/4 for the all-R term), and Hurwitz-Radon limits N anticommuting N x N operators to
  N in {1,2,4,8}; with a module of dimension M the number of operators is O(log M).  Dead end [known].
- One fixed N-dimensional space in both orders (Weyl operators of an extraspecial group, Weil representation):
  sum_{s<=k} Y_v^s X_x^s = [x=v] Id_N forces the images of the stacked X_x to be independent, so k >= N pieces [P]
  (this is the rank count with r_s <= N).  Weyl commutation X Y = omega(x,v) Y X transports a product to the other
  order only up to the phase omega(x,v), a matrix in (x,v) of rank up to N: column rank N^2 in general.
  A graded structure (dimensions from 1 to N^2) is unavoidable, by S1.  No gain found.  [P for k >= N; rest H]
- Lattice vertex operators a_x = e^{alpha_x}, b_v = e^{-alpha_v} with |alpha_x|^2 = 2h: locality with H = 2h holds and
  the top coefficient is [x=v] by Cauchy-Schwarz, so they are pieces of the SAME framework ((V) holds, growth
  Cardy): rho >= 1/sqrt2 by G9a.  [E, not computed]
- Way around (V) inside vertex algebras [O]: (V) follows from sum_s Q_s being the zero mode of a_(H-1) b in the
  vacuum line (weight wt a + wt b - H = 0).  One could ask for a_x (j) b_v = [x=v] c with c a state of positive weight
  whose zero mode is invertible on the label space (L_0 is one: scalar E0 on W_{E0}); then H = wt a + wt b - wt c.
  With c = omega this saves 2 out of 2h.  A family with wt c comparable to 2h is not known to me.

### (b) non-concave profiles
Not constructed, concavity not proved forced.  What is proved: with A4-A6 and L only nondecreasing, the staircase
profile above gives accounting cost exactly 2 lambda per piece (rho = 1/2), see section 4.  What it would need: no
states of weight in (0,h1), exactly ~N states per weight on [h1,3h1], N^2 above 3h1: two-particle states pushed from
2h1 up to 3h1.  In a vertex algebra with translation T injective on positive weight the dimensions are nondecreasing
(compatible), and lattice-type OPEs e^a(z)e^b(w) ~ (z-w)^{(a,b)} e^{a+b} do push the product state up to weight
h_a + h_b + (a,b); but the Heisenberg descendants of each e^a refill the gap.  UNKNOWN.  [H]

### (c) decision over a small field: lambda = 2 over GF(2)
Finite-N proxy used here (the right one for rho, since elementary tensors give N^2 pieces of cost 0):
  pi_k(N) = min over k-piece decompositions of max_i rank_row(A_i) * rank_col(A_i)   (rho = 1/2 means pi_k(N) = N^{2+o(1)} with k = poly(log N)).
lambda = 1: pi_2(2) = 4 = N^2 (one-bit identity, any field).  lambda = 2, GF(2), sat_k.py (exact, SAT; a SAT answer
is re-verified with GF(2) ranks, an UNSAT answer has no proof log):
  k = 4: pi_4(4) <= 16 = N^2: tensor square of the one-bit identity, profile (16,1),(4,4),(4,4),(1,16)  [N exact, tsq.py].
  k = 3: NO decomposition with all rank products <= 16.  Proof by cases (GF(2) only):
     - every profile passing the rank count is either {(<=8,<=2),(<=8,<=2),(1,*)} (or its mirror) or contains a
       column-rank-1 piece AND a row-rank-1 piece plus a middle piece (a,b), ab <= 16  [enumeration of 1100 profiles];
     - first class: UNSAT (CaDiCaL, 115 s; no proof log);
     - second class: exhaustive over all F, G (fg_exh.py, exact): min rank_row(X - f(x,u)g(y,v)) = 6 = N(N-1)/2,
       attained only by F = G = Id up to equivalence.  So the middle piece needs a >= 5 and, by the y<->u symmetry,
       b >= 5: ab >= 25.  (Consistent: SAT gives UNSAT for the middles (2,8), (3,5); (4,4) timed out in both solvers.)
     Upper bound: grouped tensor square (16,1),(6,6),(1,16), product 36 over GF(2) ((7,7), 49 over Q) [N exact].
     Hence 16 < pi_3(4) <= 36 over GF(2): already at lambda = 2, going from 4 = 2^lambda pieces to 3 pieces leaves the floor.
  k = 2: (col-rank 1)+(row-rank 1) UNSAT; a piece of column rank 1 forces the other to row rank >= 6, column rank >= 15.
  By-product: settles G9b's open case "rows = 5, cols = 1" (impossible; rows >= 6 are needed next to one column term).
Meaning: none for the asymptotic rho (the solvers could not even re-find the known 4-piece solution in 330 s, so their
timeouts carry no information); it shows exactly where the cost of grouping starts.  No GF(2) pieces with rho < 0.811
were found; char 2 keeps the status of G9b.
The full exhaustive search at lambda = 2 (2^256 tensors) is out of reach; the reduced search that is feasible is this
one: fix k and a rank profile, encode both factorisations of each piece (r_i + c_i rank-one terms, 32 variables each)
and the 256 + 256(k_both) parity constraints.

### (d) level-dependent pieces, more pieces, D >= 3 directions
- Level-dependent pieces: INESSENTIAL [P].  e_{l+1} = (rho_l + 3 zeta_l) e_l + 3/K_l, so e_l >= prod rho_j / 2 and the
  exponent is governed by the geometric mean of the rho_l; every rho_l obeys the same floor.
- More than poly(lambda) pieces: ESSENTIAL for the floor, useless for the algorithm.  Explicit example [N, exact]:
  the t-th tensor power of the one-bit identity, NOT grouped: 2^t = N pieces, each with rank_row * rank_col = N^2
  exactly (rho = 1/2); checked for t = 2 (profile (16,1),(4,4),(4,4),(1,16)).  But in Theorem A the number of terms
  |I|^{#rectangles} multiplies the time at full price, while the bits of a piece are discounted by e_l -> 0: one needs
  log|I| << e_l lambda, i.e. essentially poly(lambda) pieces.  [H: read off the cost shape of Thm A(b), not re-derived]
- D >= 3 directions: see section 3.

## 3. D directions (coordinator's question: is D^{-1/2} forced?)

Object (G6b): the D-box, 2^D corners, 2^{D-1} antipodal diagonals of lambda bits, D flattenings, each of rank
N^{2^{D-1}} =: 2^B.  rho(D) = max_piece sum_i c_i / (D B).  Floor 1/D (G6b/G6r, under their two conditions).

Accounting theorem [P, in the accounting; the hypotheses are the D-analogues of A4, A5]:
 suppose the pieces are indexed by weight vectors x in R_+^D with sum_i x_i = M, piece x is charged L(x_i) in
 direction i, the balanced point x_i = M/D is a charged piece, and |I| = 2^{o(B)}.  Then L(M) >= B - o(B) and
     rho(D) >= L(M/D)/L(M)   (= D^{-alpha} for L = c E^alpha;  >= 1/D for every concave L with L(0) = 0).
 Proof: rank count in one direction gives a piece with L(x_i) >= B - o; L nondecreasing and x_i <= M; the balanced
 piece costs D L(M/D).
So D^{-1/2} is forced exactly when alpha = 1/2 is: it is the signature of partition (Cardy, 2-dimensional chiral)
growth, not of the rank count, which allows D^{-alpha} down to 1/D at alpha = 1.
What pins alpha?  For D = 2 it is (V): M >= 4 h1 = (number of corners) x (weight of a label space) gives the second
bound 2^{alpha-1}.  The same heuristic for D directions, (V)_D: M >= 2^D h1 with L(h1) = lambda [H, my extrapolation,
NOT derived from any structure], gives  rho(D) >= max(D^{-alpha}, (2^D/D)^alpha / 2^{D-1}), minimised at
alpha* = (D-1)/D:   rho*(D) = D^{-(D-1)/D}.
   D     floor 1/D    law D^{-1/2}    rho*(D)     alpha*
   2     0.5000       0.7071          0.7071      1/2
   3     0.3333       0.5774          0.4807      2/3
   4     0.2500       0.5000          0.3536      3/4
   8     0.1250       0.3536          0.1621      7/8
Reading [H]: the accounting does NOT force D^{-1/2} for D >= 3; it would be beaten by a graded space with
log dim ~ E^{(D-1)/D} (the growth of plane partitions for D = 3, of a free field in D space-time dimensions in
general) carrying a D-direction locality identity.  No such identity is known to me; with partition growth
(alpha = 1/2) the law D^{-1/2} is exactly what the accounting gives.  This is a target, not a result.

Smallest case, exact [N]: lambda = 1 (N = 2), class of "epsilon-matching" terms.  Up to an invertible change of basis
at half of the corners (ranks unchanged) the D-box is T_eps = prod over antipodal pairs of eps(w_p, w_pbar).  A term
= product of eps over a perfect matching of the 2^D corners; its rank in direction i is 2^(number of pairs cut by
coordinate i), so its total cost is the sum of the Hamming distances of its pairs (box: D 2^{D-1}; floor 2^{D-1},
attained exactly by matchings of the cube GRAPH Q_D by edges).
   D = 2: T_eps = difference of the 2 edge matchings (one-bit identity): max cost 2 of 4, "rho" = 1/2 = floor.
   D = 3: T_eps is NOT in the span of the 9 edge matchings of Q_3 (rank 9 vs 10, exact over Q);
          it IS a +-1 combination of 10 matchings of cost 4 or 6 (exact identity, sympy): max cost 6 of 12,
          "rho_3" = 1/2 in this class; 1/3 < 1/2 < 3^{-1/2} = 0.577.
   D = 4: T_eps is NOT in the span of the 272 edge matchings of Q_4 (rank 272 vs 273, exact, flint);
          it IS in the span of 16 matchings of cost 16 (two commuting one-bit expansions): "rho_4" <= 1/2 = 4^{-1/2};
          NOT in the span of all 10064 matchings of cost <= 10 [Np: rank 1262 vs 1263 mod 2^31-1 on 3000 coordinates];
          the matchings of cost <= 12 span the whole invariant space (rank 1430 = Catalan(8) mod p, hence over Q),
          so T_eps is in their span: min max cost = 12 of 32, "rho_4" = 3/8 in this class; 1/4 < 0.375 < 4^{-1/2} = 0.5.
Caveats: (i) these are numbers of the UNGROUPED smallest case (number of pieces not constrained), the same sense in
which D = 2 has "rho = 1/2" at lambda = 1 while the grouped tensor powers give 0.896 (Q) / 0.811 (GF(2)) and G6b's
grouped D = 3 search gives 0.914; they say nothing about asymptotic rho(D).  (ii) The class is restricted: terms need
not be SL_2-invariant; e.g. GHZ x GHZ on two opposite faces has cost (2,2,0) = floor and is not a matching.
What the smallest case does show: for D = 3, 4 the floor 1/D is not reached by the skein mechanism that reaches it
for D = 2; the direct D-box is harder than iterating D = 2, in agreement with G6b.

## 4. Computations (all under 10 minutes each)
 file            what                                                              result
 profiles.py     accounting floors of section 2 (max_x L(x)+L(M-x))/(4 lambda)      0.7071 / 0.5 / 1.0 / 0.5 / 0.6667 (profiles.out)
 tsq.py          tensor square of the one-bit identity, N = 4, exact ranks Q, GF(2)  (16,1),(4,4),(4,4),(1,16); grouped middle (7,7) Q, (6,6) GF(2)
 sat_k.py        GF(2), N = 4, k pieces with prescribed rank profile (Glucose)      k=2 (.,1)+(1,.) UNSAT; k=3 middle (2,8) UNSAT 3 s; (4,4), (8,2)^2, k=4: timeout 240 s
 sat_kc.py       same with CaDiCaL (no interrupt; external timeout 330 s)           k=3 (8,2),(8,2),(1,.) UNSAT 115 s; middle (3,5) UNSAT 13 s; middle (4,4) and k=4 (known SAT): no answer
 fg_exh.py       exhaustive GF(2): min rank_row(X - f(x,u)g(y,v)), N = 4            rank F = 0..4: 16, 15, 13, 10, 6;  m* = 6  (fg_exh.out)
 dmatch.py       D = 2, 3 box, lambda = 1, eps-matchings, exact over Q              D=3: not in span of 9 edge matchings; 10-term identity with costs <= 6
 dmatch4.py      D = 4, edge matchings of Q_4 (272), exact over Q (flint)           not in span (272 vs 273); in span of 16 matchings of cost 16
 dmatch4c/d.py   D = 4, costs <= 10, <= 12, mod p                                   <= 10: no (all 10064, 3000 coords); <= 12: spans all 1430 invariants
Outputs: sat_k*.out, sat_c*.out, dmatch4c.out, dmatch4d.out, fg_exh.out, profiles.out.  The first dmatch.py / dmatch4.py
outputs are quoted in section 3 (rerun: about 10 s and 60 s).

## 5. Verdicts (ESSENTIAL = explicit example below the floor when dropped; INESSENTIAL = the proof extends; UNKNOWN)

 assumption                     verdict       why
 A1-A3 one graded space,        INESSENTIAL   as a restriction it is empty (Prop. 1 [P]): every piece family is of this form with H = 0.
   both-orders identity                       The floor does not come from it.
 grading vs filtration          INESSENTIAL   the proof uses only L(E) = log dim W_{<=E} [P].
 A4 dimension accounting        (tight)       for App. B the true ranks equal the window dimensions (G18c [N]); no slack to exploit there.
 A5 interval / balanced piece   UNKNOWN       needed by the proof (without it the accounting floor is 1/2: pieces only at the two ends);
                                              no family with empty middle is known.
 A6 (V): H >= 2 h1, M >= 4 h1   UNKNOWN,      without it the accounting floor is 2^{-alpha} -> 1/2 (linear L, H = o(E0)) [P, accounting];
                                proof-        NO operators found.  It is a theorem inside vertex algebras when the identity is the vacuum
                                essential     coefficient [E]; a non-vacuum zero mode would shorten H by wt c only [O].
 A7 concavity / power law       UNKNOWN,      the staircase profile gives accounting rho = 1/2 even with (V) [P, accounting]; not
                                proof-        constructed, not excluded.
                                essential
 A8 field                       INESSENTIAL   for the floor (rank count over any field); App. B needs char 0 or > h.  Fixed char 2: best
                                (floor)       known 0.811 (G9b); lambda = 2 with 3 pieces leaves the floor [N, section 2c].  rho_F < 1/sqrt2: UNKNOWN.
 A9 counting vs decision        INESSENTIAL   (G9b: decision = fixed small field = A8.)
 A10 poly(lambda) pieces        ESSENTIAL     2^lambda ungrouped tensor-power pieces have rho = 1/2 exactly [N exact at lambda = 1, 2; P in general:
                                (but forced   ranks multiply].  Useless: the number of pieces is paid at full price in Theorem A [H].
                                by Thm A)
 A11 same pieces at each level  INESSENTIAL   the exponent is the geometric mean of the rho_l [P].
 A11' D = 2 directions          UNKNOWN       accounting: rho(D) >= D^{-alpha}; D^{-1/2} is forced only if alpha = 1/2 is; the speculative
                                              optimum is D^{-(D-1)/D} [H].  Smallest case, matching class: 1/2 (D=2), 1/2 (D=3), 3/8 (D=4):
                                              the floor 1/D is reached only for D = 2 [N exact for D = 3; Np for D = 4].

Best rho found: NO piece family with asymptotic rho < 1/sqrt2 (none even below the known 0.811 over GF(2) or
1/sqrt2 over Q).  Values below 1/sqrt2 exist only (i) as accounting profiles without operators (0.5), (ii) with
exponentially many pieces (0.5, exact, useless), (iii) in the ungrouped smallest cases (exact, not asymptotic).

Consequence for the exponent.  Unchanged: 2^{O(m/sqrt(log m))} with rho = 1/sqrt2; the only proved floor for general
pieces is 1/2, i.e. the method is not shown to stop before 2^{O(m/log m)}, and cannot go below that.  What changed is
the diagnosis: "one graded space" is not the obstacle; the obstacle is the pair (V) + Cardy growth of vertex algebras.
A gain needs either (1) N x N operator families on a graded space of exponential growth that commute up to a window
H = o(E0) [would give rho -> 1/2], (2) a staircase profile with (V), or (3) a D-direction locality identity on a space
with log dim ~ E^{(D-1)/D}.  (1) and (2) would in particular improve the Hurwitz-type bound P(N) <= N^{sqrt2+o(1)}
(rho < 1/sqrt2 implies P(N) <= N^{2 rho + o(1)}), so they are at least that hard; for (3) I know no such link.

## 6. Adversarial self-review
- Prop. 1 is trivial; its only content is to locate the floor in (V) + shape.  It does not produce pieces.
- All "floors 1/2" of section 2 are statements about the ACCOUNTING (upper bounds on ranks by dimensions); they show
  where G9a's proof stops, not that better pieces exist.  No verdict ESSENTIAL is claimed for (V) or concavity.
- (V)_D and the table rho*(D) = D^{-(D-1)/D} are an extrapolation of mine (M >= 2^D h1), not derived from any algebraic
  structure; the "free field in D dimensions" remark is an analogy.  The D-direction accounting theorem assumes a
  balanced charged piece and conservation sum x_i = M, as G18c's Hermite picture suggests; unverified for real schemes.
- Smallest-case numbers (D = 3: 1/2; D = 4: 3/8; lambda = 2: pi_4 = 16 < pi_3) are ungrouped/finite and restricted to
  the eps-matching class for D >= 3; they do not bound asymptotic rho(D) in either direction.  The D = 4 statements
  "not in the span at cost <= 10" is mod one prime on 3000 of 12870 coordinates (the direction "not in span mod p" is
  rigorous mod p; over Q it is very likely but not proved).  "Cost <= 12 suffices" uses dim(invariants) = Catalan(8).
- lambda = 2, k = 3: the first class rests on one UNSAT answer of CaDiCaL without a proof log; the second class is an
  exhaustive enumeration (my own rank routine, cross-checked only through agreement with the SAT answers (2,8), (3,5)
  and with rank 6 = N(N-1)/2 at F = G = Id).  Over Q nothing was proved for k = 3.
- "More pieces are useless" was read off the shape of Theorem A(b), not re-derived.
- The (V)-is-a-theorem-in-vertex-algebras remark and the lattice-vertex-operator remark are sketches [E].
- Not done: an exact rank computation of App. B pieces in the regime where the window bites (h >= 16, N >= 231:
  53361 x 53361 matrices), so "A4 is tight" rests on G18c's small cases h <= 5.
- Nothing here is a lower bound for 3-SAT or an improvement of any algorithm.
