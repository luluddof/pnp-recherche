# G18c — Hybrid of the two October-2026 sub-exponential 3-SAT algorithms (working notes, unreviewed)

Status tags: [P] proved here (short argument written), [N] numerical, [H] heuristic, [O] open.
(A) = Claude/Anthropic recursion (pieces of the crossing; App. A + App. B).
(B) = OpenAI App. D (Cauchy deformation, tapered Hermite, 2-adic modes).
Neither paper's correctness is assessed here; everything is "if (A), (B) are right as written".
Sources read: claude_3sat.txt (sec. 2, App. A, Def. 7, Thm A), OpenAI paper.txt App. D.4, our notes
G9a, G9b, G14b, G17a, G17c.

## 0. Intermediate save (written before the numerics; refined below)

### Common formalism (draft)
Object: a closed tensor network = vertices with arbitrary tensors + M one-bit wires; a wire is the
identity [Y=X] between its two ends. Fix a LABELLING of the vertices by words c = (c_1..c_D) over an
alphabet; coordinate i defines the partition pi_i ("direction i") of the vertices.
One level = an identity   c0 * (identity on the wires) = sum_k kappa_k tau_k,
each tau_k "cheap in one direction i(k)": after merging the classes of pi_i(k), tau_k has few bits
between classes. Both algorithms are HERMITE INTERPOLATIONS of a one-parameter deformation of the
identity kernel between direction-local nodes, and the jets are priced by a partition (Fock) count.

## 0b. Numerics done so far (saved before the write-up)  [N]
- out_A_profile.txt (heisS_copy.py = our G9a script + one print): Claude App. B pieces, h = 3,4,5
  (N = 3,5,7).  rank_row(A_s) as a function of x = s+2h: 1,2,4,7,12,18,27,40,48,46,49,... i.e.
  = sum_{w<=x} p(w) (1,2,4,7,12,19,30,...) until it saturates at N^2; rank_col mirrored; A_{-2h} has
  row rank 1, A_{+2h} has column rank 1.  => App. B is a TWO-NODE Hermite scheme (nodes T=0: rows,
  T=infinity: columns; target T=1) for J(T) = sum_s T^{s+2h} A_s, jet x priced by dim W_{<=x}.
- out_B_a1.txt, out_B_a2.txt (bkernel.py): OpenAI kernel K(T) on one packet, L = 2, 4, exact mod a
  61-bit prime.  K(0) = I; identity has rank 2^L in every direction (L=4 packet = Claude rectangle,
  N = 4); K(t) has rank 1 in direction v2(t); BUT the exact jets [z^j]K(t+z) have rank 6 at j = 1 and
  FULL rank 16 for 2 <= j <= 12.  Over Q, OpenAI's jets are not cheap pieces.
- out_B_smith.txt (bsmith.py): Smith form over Z_(2) of the flattened K(t + alpha*zeta).  Rank over Q
  is full (16) but the 2-adic valuations are graded:  dir 1 (t = +-2): {0:1, 2:4, 4:6, 6:4, 8:1}
  (= (1+q^2)^4, weight 2 per crossing bit);  dir 0 (t odd): {0:1,3:2,5:2,6:1,7:2,8:2,10:3,12:2,14:1}.
  So "rank mod 2^p" N_p = 1,3,5,6,8,10,13,15,16 grows with the precision p: the 2-adic weight
  filtration is OpenAI's counterpart of Claude's weight window.

## 1. Both algorithms in one formalism

Object. Closed tensor network: vertices with arbitrary tensors, M one-bit wires, wire = identity
[Y=X] between its two ends.  (A)'s (K,lambda)-clique is this with K vertices and lambda parallel
wires per pair; (B)'s Boolean bipartite network is this directly.
Labels. Vertices get words c = (c_1..c_D) over an alphabet; coordinate i gives a partition pi_i of the
vertices ("direction i").  A GADGET is a set of wires that is a full crossing in every direction:
(A) rectangle = two diagonals of a 2x2 sub-grid (D = 2: rows, columns); (B) packet = L = 2^a slots, in
every coordinate the unordered symbol pair is fixed and bit i of the slot gives the orientation.
[N] a packet with a = 2 IS a rectangle with lambda = 2 (identity rank 16 = N^2 in both directions).
One level. A deformation J(T) of the gadget's identity kernel, polynomial in T, equal to the identity
at the target point and direction-local (rank 1 across the classes of pi_i) at the nodes; the identity
is recovered by HERMITE interpolation; the jet of order x at a node is not local any more and its
cost is the dimension of a weight-truncated graded (Fock-type) space.

|                         | (A) Claude, Thm A + App. B pieces            | (B) OpenAI App. D                                  |
|-------------------------|-----------------------------------------------|-----------------------------------------------------|
| deformation             | J(T) = sum_s T^{s+2h} A_s, degree 4h (*)      | K(T) = prod_x l_{Y_x}((W(T)X)_x), W Cauchy, deg L(2L-1) per packet |
| variable                | one T_r PER rectangle (product of sums)       | ONE global T (degree M(2L-1))                       |
| target                  | T = 1: sum_s A_s = h! * crossing              | T = 0: W(0) = I                                     |
| nodes -> direction      | T=0 -> rows, T=infinity -> columns (2 nodes)  | t in {-L..L-1}\{0} -> coordinate v2(t), u <= a used |
| exactness at node       | A_{-2h} row-rank 1, A_{+2h} col-rank 1 [N]    | W(t) = c_t * shift, rank 1 [N]                      |
| jets used               | x = s+2h <= 2h at each node (exactly determined) | q_t = h*(u - v2(t)), h = 3M/u (tapered, ~3x over-determined) |
| what prices a jet       | weight window: rank <= dim W_{<=x} = sum_{w<=x} p(w) [N] | 2-adic precision p_t ~ M D_t/u: tuples with sum w_Z beta_Z < p_t |
| the graded space        | Q[t_1,t_2,..], one oscillator per weight      | modes Z=(packet,tile,b), weight b+1, 6M/D_t per weight |
| "wire" between classes  | theta in a basis of W_E (c_row bits)          | mode tuple beta                                     |
| what is done with it    | KEPT as wires -> smaller clique -> recursion  | ENUMERATED (exp(O(M/sqrt u)) terms), children disconnected |
| # terms per level       | |I|^{|R|} = (2 lambda)^{O(K^2)} (negligible)  | 2^{O(M/sqrt u)} (this IS the cost)                  |
| level = on the labels   | K = g^2 -> g: number of coordinates HALVED    | all u directions at once; children: M/8, u/2        |
| where the time goes     | leaf: brute force on rho^l S bits (+ S/K sums)| root: number of tuples                              |
| per-level gain          | rho = 1/sqrt 2 with 2 directions              | "rho_B" = c/sqrt u with u directions                |
| analogue of the other's | taper <-> window split 2h | 2h between the two nodes; precision <-> max weight 4h | rho <-> log2(#tuples)/M |
| label resource          | D = log K / log K_0, K <= m^{1/9}             | u <= a, r^{2a} <= M/(2 Delta)                        |
| result                  | e = Theta(D^{-1/2}) = Theta(1/sqrt(log m))    | Theta(u^{-1/2}) = Theta(1/sqrt(log s))               |
| ring                    | Z (or GF(p), p > h); divides by h!^{|R|}      | Z_(2)-order mod 2^p; divisions = 2-adic losses      |
(*) reading of App. B through our G9a axioms: A_s factors row-wise through weights <= s+2h, column-wise
through weights <= 2h-s.

Where each gets its sub-exponential saving, in the same words:
 the 2^{bits} states of the crossing wires are embedded in a graded space whose dimension up to weight E
 is exp(Theta(sqrt(E * species))) (Hardy-Ramanujan in (A) - "logarithms of sums of partition numbers";
 the generating product prod (1-e^{-theta j})^{-6M/D} of (D.15) in (B)).  Full resolution needs the whole
 weight range R (A: 4h for 2 lambda bits with lambda ~ sqrt h; B: ~ M*D for M bits, i.e. D^2 per tile of
 D bits).  Interpolation lets each DIRECTION carry only its share of the jets, i.e. weight <= R/D_dir
 (A: 2h = R/2, two nodes; B: p_t ~ M D/u = R/u, the taper).  With sqrt-growth the cost is
 bits * sqrt(1/D_dir):   rho(D_dir) = D_dir^{-1/2}.     (LAW)   [H; exact for (A) at D_dir = 2; O(.) for (B)]
 (A) reaches many directions hierarchically (2 per level, l levels, D = 2^l coordinates consumed),
 (B) flat (u directions in one level).  In both D <= O(log M) because every gadget needs populated
 pair classes (2^D orientation classes of wires per gadget; K^2 <= m^{2/9}; r^{4a} <= M/Delta).

## 2. Substitutions

Composition rule used below [H]: a level that uses D_1 coordinates leaves labels with D/D_1
coordinates, so a tree of levels with direction counts D_1, D_2, .. has prod D_j <= D = O(log M)
and gain prod D_j^{-1/2} >= D^{-1/2}.  The savings MULTIPLY, the resources multiply too, and x -> x^{-1/2}
is multiplicative: every composition lands on the same exponent c*M/sqrt(log M).

(i) Claude's pieces as kernels/wires of OpenAI's contraction.
 (i-a) pieces instead of the Cauchy kernel, per packet: well-defined for a = 2 packets; it IS
   algorithm (A) (packet = rectangle).  NO GAIN (identity).
 (i-a') same pieces with ONE global T (OpenAI-style): F(T) = sum_sigma T^sigma Z_sigma, nodes 0, infinity;
   the row-merged jet x is priced by an n-coloured Fock space (n = rectangles between two rows) at
   weight x; balanced x = 2hn: log = pi sqrt(2 n x /3) = n * pi sqrt(2(2h)/3), the SAME as n independent
   windows [P, from the Meinardus/Hardy-Ramanujan leading term].  Global budget = no change.  NO GAIN.
 (i-b) Claude's recursion on OpenAI's beta ("modes as wires", listed as open in G17a): beta is
   a wire between the two classes of a packet, but the saving exp(sqrt(n p)) exists only for the
   JOINT budget sum_Z w_Z beta_Z < p_t.  As independent per-packet wires the alphabets have
   (3M/L) * sqrt(2 L M/u) = 3M sqrt(2M/(L u)) >> M bits (L <= M^{1/24} by (D.6)); Claude's recursion gains a
   factor e >= 1/sqrt(log) on bits, never sqrt(L/M).  Its balancing/padding step also forgets the
   grading.  ILL-DEFINED as a gain (well-defined as an algorithm, strictly worse).
(ii) OpenAI's deformation + interpolation inside Claude's level (instead of the piece identity).
 Over Q, per gadget: well-defined, pieces = jets [z^j]K(t+z) with an exact Hermite identity, but
 [N] at L = 4 the jets have rank 6 (j = 1) and FULL rank 16 (j >= 2) in their own direction, and a
 per-packet identity needs about 5 jets per node: rho = 1.  NO GAIN.  Reason: (B)'s jets are cheap
 only modulo 2^p (Smith valuations graded, [N]); over a field they are not pieces in the sense of
 Definition 7 / our "F-pieces".
 Mod 2^p with a global T: this is (B) itself restricted to the 2 directions of a grid level:
 "rho_B(2)" = c/sqrt 2 with c > 1 from (D.15) (my estimate: c in the tens with the paper's constants, since p_t >= 24 M D_t/u and there are 6M/D_t modes per weight), against exactly
 1/sqrt 2 for App. B.  As a u-direction level inside (A)'s recursion (K = g^u -> g): gain u^{-1/2},
 coordinates divided by u: by the composition rule, NEUTRAL.  NO GAIN in the exponent.
(iii) 2-adic truncation applied to Claude's recursion.
 Well-defined (the recursion is ring-agnostic, G9b; d = h! must be handled by precision > v2(d^{|R|})).
 [P] It cannot lower rho below 1/2, nor remove the window rank count of G9a: if sum_i A_i =
 d*X + 2^p E with p > v2(d), then det(d X + 2^p E) != 0 over Q, so sum_i rank(A_i) >= N^2 in both
 flattenings whatever the truncation (same argument as G14b Thm A); and a rank-one decomposition
 mod 2^p lifts to one over Z of a matrix congruent to d X.  What truncation buys in (B) is the price of
 a Taylor order Theta(M) for a GLOBAL variable; (A) has no global variable and no Taylor order to
 price (its per-rectangle expansion is finite and exact).  NO GAIN.
(iv) One algorithm as base case / sub-call of the other.
 (iv-a) (B) as base case of (A) at level j*: needs (B) with Delta = 2M/K instead of s^{9/10}, hence
   r^{4a} <= K_{j*}/4, u <= log K_{j*}/48 (UNPROVED extension of Thm D.1; entry cost 2^{2N/K} must stay
   below (B)'s allowance).  Exponent: rho^{l-j*} * c_B/sqrt(u) = c' / sqrt(2^{l-j*} log K_{j*}) =
   c'/sqrt(log K_l), INDEPENDENT of j*.  NEUTRAL: no gain, only the constant can move.
 (iv-b) (A) at the leaves of (B): (B) is root-dominated (exponents sum geometrically,
   sum_d M_d/sqrt(u_d), the leaf 2^H is tuned to the root term), so a cheaper leaf changes nothing;
   stopping (B) at depth d and calling (A) gives sum_{j<d} M_j/sqrt(u_j) + M_d e_A >= the root term.
   NO GAIN.

## 3. Exponent of the most promising hybrid (heuristic)
 Hybrid H = (A) on the top l-j* levels, a (B)-type flat level with u directions below (or (iv-a)).
 exponent(H) = c * m / sqrt(prod of direction counts) with prod <= C log m, i.e.
     2^{O(m / sqrt(log m))}   (sparse case; n loglog n / sqrt(log n) in general) -- THE SAME as (A) and (B).
 Status: heuristic neutrality (LAW + composition rule); nothing better than either paper is claimed,
 and no constant was computed.  With rho(D) = D^{-gamma} the same computation gives m/(log m)^gamma, so
 a gain needs a level with rho(D) < D^{-1/2}: e.g. the floor rho = 1/2 at D = 2 would give m/log m.
 Key identities tested numerically: section 0b (three files).

## 4. Verdicts
 (i-a)  NO GAIN  - it is (A).                       (i-a') NO GAIN - global budget = product of windows.
 (i-b)  ILL-DEFINED as a gain - beta is budget-coupled, not a product of pairwise wires.
 (ii)   NO GAIN  - over Q the jets are full rank [N]; mod 2^p it is (B) with D = 2, same 1/sqrt 2 law.
 (iii)  NO GAIN  - rank floors survive 2-adic approximation [P]; (A) has no Taylor order to price.
 (iv-a) NO GAIN (neutral for every switching level; needs an unproved extension of Thm D.1).
 (iv-b) NO GAIN  - (B) is root-dominated.
 COMMON REASON: both algorithms are the same phenomenon - Hermite interpolation of the identity
 between direction-local kernels, with jets priced by a partition count; both obey
 exponent = M / sqrt(#label coordinates), #coordinates = O(log M).  Their savings draw on the SAME
 resource (label coordinates / populated pair classes) through the SAME growth law (sqrt), so they
 multiply only by splitting the resource.  They differ in bookkeeping: per-gadget variable + recursion
 on the residual wires (A) versus global variable + 2-adic pricing + enumeration (B).

## 5. What would be a real gain (open)  [O]
 - a graded space / kernel family with log dim ~ E^alpha, alpha < 1/2, compatible with an interpolation
   identity (G9a excludes it in the vertex-algebra window geometry; not excluded in (B)'s 2-adic
   geometry: one tile has D bits and weight range ~D^2 - is that forced?);
 - D-direction pieces for (A) (a D-point locality identity on a D-cube gadget) with the per-term
   "min" accounting of (B): floor 1/D, would give m/log m if attained;
 - a proof of rho(D) >= c D^{-1/2} for the whole family (G14b's direction-rank problem) would turn the
   neutrality into "sqrt(log) is optimal for every hybrid".

## 6. Adversarial self-review
 - The "two-node Hermite" reading of App. B is partly a tautology (any graded family of pieces is the
   coefficient list of a Laurent polynomial).  The non-trivial content is the rank profile
   sum_{w<=x} p(w) and rank 1 at both ends, checked only for h <= 5 (N <= 7), where the pieces are far
   from efficient (rho_eff = 1 at these sizes).
 - LAW for (B) is an O(.) statement with a constant > 1 at small u; "rho_B(2) vs 1/sqrt 2" compares an
   asymptotic law with an exact value.  The R/D_dir "share of the weight range" explanation is my
   paraphrase of (D.10)-(D.15), not a theorem.
 - The composition rule (resources multiply) is exact for (A)'s grid levels; for a flat level followed
   by further levels it assumes r = r'^{a'} labels and was not built.
 - (iv-a) relies on running Thm D.1 with max degree M/K^{O(1)}; not checked.  (i-b) compares leading
   terms only.
 - Numerics: one packet, L <= 4 (a = 3 attempted, see below), ranks mod a 61-bit prime (lower bounds
   for Q-ranks; they came out full, so exact), Smith forms at zeta = 1, 3 only; no joint-budget test
   with several packets.  bkernel's rank-1 check at t = +-3 is rank 1 in BOTH directions (only one
   slot pair survives the shift), which is a small-L artefact.
 - Not explored: hybrids that are not substitutions (multi-point Heisenberg kernels, odd-prime grids,
   Claude pieces over a 2-adic Fock space with weighted bases).  "Same phenomenon" is a structural
   claim supported by the dictionary and tiny numerics, not a theorem.
 - Neither paper's correctness was assessed.  Nothing here is a lower bound for 3-SAT and nothing
   improves either published exponent.

## 7. Late addenda
 - a = 3 (L = 8) Smith profile: attempted with a 280 s limit and TIMED OUT (out_B_smith_a3.txt contains
   only "exit 124"): there is NO a = 3 data; all (B) numerics are for one packet with L = 2 and L = 4.
 - Precise form of the profile check on the column side (out_A_profile.txt, h = 5):
   rank_col = 1,2,4,7,12,19,30,39,48,49 for 2h-s = 0..9, equal to sum_{w<=x} p(w) up to x = 6;
   the row side is slightly smaller from x = 5 on (18, 27, 40): "<= dim W_{<=x}", not "=".
 - All background processes were killed before reporting.
