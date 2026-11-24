# G19a — the "Hamming-distance-1 loads" opening toward 2^{O(n/log n)} (model DF, App. D of the OpenAI preprint)
(AI agent notes; interim save 0.  Levels: [P] proof written here, not reviewed; [N] numerical; [C] conjecture; [E] evidence.)
Sources read: g18s/notes.md, g18e/notes.md, g18d/notes.md, g14b/notes.md; preprint paper.txt lines 3007-3060, 3290-3362, 3560-3600.
Plan: 1 construction + exact check; 2 origin of the load constraints; 3 separator inequality + expander bound;
4 intermediate formula; 5 verdict.

## [save 1] 1. The non-uniform construction, exactly   [P + exact check c1_construct.py / c1_out.txt]
Setting (DF): network = multigraph G, M wires, wire e has two end bits (Y_e, X_e) on its two end vertices;
I = prod_e [Y_e = X_e].  Directions i = 1..u = partitions pi_i of the vertices (App. D: group = symbol of the label in
coordinate i).  A wire is CUT by direction i if its ends lie in different groups of pi_i;  M_i := #wires cut by i;
d_e := #directions cutting e (= Hamming distance of the two end labels on the u direction coordinates).
"Loads on Hamming-distance-1 pair classes" = every wire has d_e <= 1 (its two end labels differ in at most one direction
coordinate).  Then sum_i M_i <= M, hence min_i M_i <= M/u.
CUT-CONDITIONING IDENTITY (valid for EVERY labelling, not only Hamming-1):  for any direction i
     I = sum_{z in {0,1}^{cut wires of i}}  prod_{g in pi_i} tau_g^z ,
     tau_g^z = prod_{e uncut, inside g}[Y_e = X_e] * prod_{e cut, end bit b of e located in g}[b = z_e].
 Exact over Z (no 2-adic error), kappa = 1, each term is a product over the groups of ONE direction: a DF identity with
 N = 2^{M_i} terms.  It is nothing but conditioning on the wires of an edge cut (the classical cut-set / separator step).
 So in DF, for every labelling,     log2 N  <=  min_i M_i.                                              (UB)
Lower bound: A'' with ALL cut wires of ... : if every quotient multigraph G/pi_i is bipartite (always true for r = 2) and
 d_e <= 1, orient each cut wire from the "source" side to the "sink" side of its unique direction; then for direction i
 (1/2) sum_g |e1_i(g) - e2_i(g)| = M_i (uncut wires cancel inside their group), so log2 N >= min_i M_i:  N = 2^{min_i M_i}.
 For r >= 3 with non-bipartite quotients the trivial identity is NOT optimal: one direction, 3 groups, one wire per pair
 of groups is the 2x2 matrix-multiplication tensor <2,2,2>, rank 7 < 8 (Strassen; recalled, not re-derived); A'' gives 4.
Exact check (own code, integers, full tensors on 2M bits, M <= 6): for 6 instances the 2^{M_i}-term sum equals I for
 every direction, every tau_g is constant in the bits outside g, and brute force over all 3^M oriented subsets gives
 A'' = min_i M_i on the four Hamming-1 instances (u=2: M_i=(2,2),(3,2); u=3: (2,2,2),(1,2,2)), with rank certificates
 (rank of restricted I = 2^w, rank of restricted terms <= 2^{w-A''}).  263 random labellings (u = 2,3, M <= 6):
 A'' <= min_i M_i always (0 violations), equality in 215.
GENERAL SANDWICH [P]:  A''(labelling) <= log2 N_DF <= min_i M_i  for every labelling.  The "2^{M/u}" of G18s 2c is
 the upper bound (UB) evaluated on a labelling with sum_i M_i <= M.  It is an upper bound of the SEPARATOR type.

## 2. Where the load constraints come from   [reading of App. D.4.2, p.58-59]
Combinatorial statement.  Network = bipartite multigraph (rows, columns), M wires, max degree <= Delta = ceil(s^{9/10}).
Labelling = map lab: V -> [r]^a.  Pair class gamma = (row label, column label) in [r]^a x [r]^a (D_arr = r^{2a}
classes); load X_gamma = #wires whose row end has the first label and column end the second.
 * The labelling is chosen by the DESIGN: "Independent uniform vertex labels give E sum_gamma (X_gamma - M/D_arr)^2 <=
   2 M Delta ... Conditional expectation gives a labeling with this bound" (Lemma D.3; derandomised hashing).  Nothing in
   the instance forces it: the lemma holds for EVERY network with M >= 2 D_arr^2 Delta, M >= L D_arr, 2L <= Delta (D.6).
 * (D.6) is a PARAMETER inequality (how large a, hence u <= a/2, may be for given M and Delta); the instance enters only
   through M and the degree bound (which is why variables get equality trees: bounded degree makes hashing uniform).
 * Why the design WANTS uniform loads: packets must contain, for a fixed unordered symbol pair in every coordinate,
   all 2^{a'} orientations equally often ("using 2^{a-a'} slots from each of these labels per packet"), so that
   slot-number shift by t is local in coordinate v2(t) (Cauchy kernel).  Uniformity is needed by the Cauchy/packet
   machinery, not by the model DF.  A designer who gives up uniformity also gives up App. D's identity: outside
   near-uniform loads the only identity available today is (UB).
So the question "can the designer choose other loads on a hard instance" is legitimate; answer in sec. 3.

## [save 2] interim findings (details and proofs in sec. 3-5 below)
 * c2_worst.py (search for the worst label-class measure): u=2: V=1/2; u=3: V=1/3 found (three of the four antipodal
   classes).  Led to the STAIRCASE / CIRCULAR labelling (2u classes on a cycle, directions = the u diameter cuts) for
   which the flattening bound A'' is <= M/u on EVERY network (discrete intermediate-value argument), although every
   direction is a balanced cut.  So A'' cannot prove M/sqrt(u) for arbitrary labellings.
 * Lemma Q (quadratic local optimum) gives the matching general lower bound A'' >= min_i M_i/u - u (r = 2).

## [save 3] 3. The key question: can a hard (expander-like) network be labelled with concentrated loads?
Notation: vol(S) = number of wire ends on S (vol(V) = 2M); cut(S) = #wires leaving S; conductance
phi = min_{vol(S) <= M} cut(S)/vol(S).  A direction is beta-BALANCED if every group has vol <= beta*2M (App. D:
Lemma D.3(3), beta = 3/r; needed so that children are smaller).  phi' := phi * min(1, 2(1-beta)).

(3.1) Counting identity [P].  sum_i M_i = sum_e d_e.   (each wire is counted once per direction cutting it)
(3.2) Expander cut bound [P].  On a network of conductance phi every beta-balanced direction has M_i >= phi' M.
   Proof: 2 M_i = sum_g cut(g).  Groups with vol <= M give cut >= phi vol.  If one group g0 has vol > M then
   cut(g0) >= phi (2M - vol g0) >= 2 phi (1-beta) M and the other groups together also give >= phi (2M - vol g0).  QED
(3.3) Separator inequality [P].  If a fraction >= 1-eps of the wires have d_e <= 1 (labels differ in at most one of the
   u direction coordinates), then some direction has M_i <= M (1/u + eps (1 - 1/u)); if that direction is beta-balanced,
   merging its groups gives a bipartition with both sides of volume >= min(1-beta, 1/3) * 2M cutting <= M(1/u + eps)
   wires (a balanced EDGE separator; a vertex separator of at most that size).  Proof: (3.1) with d_e <= u.  QED
   Equivalent form: "all wires at distance <= 1" = the label map is a homomorphism of the network into the reflexive
   Hamming graph H(u,r).  Edge-isoperimetry of H(u,r) itself is irrelevant (classes may carry unequal loads); what is
   forced is only the averaging (3.1).  The same inequality holds for any distance profile: min_i M_i <= M * mean(d_e)/u.
(3.4) Consequence on hard instances [P].  Conductance phi, balanced directions:
     mean d_e >= phi' u;   fraction of wires with d_e <= 1  <=  1 - (phi' u - 1)/(u - 1)   (-> 1 - phi');
     Hamming-1 loads (eps = 0) are possible only if u <= 1/phi';
     the cut-conditioning identity (the ONLY identity behind "2^{M/u}") has 2^{min M_i} >= 2^{phi' M} terms for every
     labelling.  The opening of G18s 2c is therefore closed on expanders: u cannot grow, and the count is 2^{Omega(M)}.
   Numbers: random 3-regular graphs have bisection width >= n/9.9 (Kostochka-Melnikov; recalled), i.e. >= 0.067 M;
   c3_expander.py: 4-regular n = 48: Cheeger bound phi >= 0.086, observed min_i M_i/M in [0.40, 0.54] for random balanced
   labellings.  NOT computed: phi of the group graph of a random 3-CNF after App. D's subdivision/equality trees
   (subdividing edges and bounded-degree equality trees keep conductance Theta(1), constant unverified).
(3.5) What else is paid.  In DF the classes are not tensors; the children are the groups (sub-networks).  A recursion
   that uses only cut-conditioning identities is a recursive edge-separator (branch-decomposition) contraction: cost
   >= 2^{cut at the root} >= 2^{phi' M}, i.e. 2^{Omega(n)} for 3-SAT at constant density (M = Theta(n)); if instead the
   directions are unbalanced (one group with volume (1-delta) 2M), the cut can be small (>= 2 phi delta M) but that
   child keeps >= (1-2 delta) M wires: no progress.  So yes: large classes/groups force 2^{Omega(n)} in this scheme.

(3.6) General lower bound for ARBITRARY labellings [P]  (new; "nothing is proved" -> M/(2u)).
   Lemma Q.  r = 2.  For every labelling there is an oriented wire subset with  min_i b_i >= min_i M_i/u - u,
   b_i = (1/2) sum_g |e1_i(g) - e2_i(g)|.  Hence  log2 N >= min_i M_i/u - u.
   Proof.  Fix a reference orientation of each wire and delta_e in {-1,0,1}^u, delta_e(i) = +-1 iff direction i cuts e
   (sign = side of the head).  For x in [-1,1]^E, b_i(x) = |sum_e x_e delta_e(i)| (x_e = +-1 selects e with an
   orientation).  Let s in {+-1}^u maximise Q(s) = sum_e <s, delta_e>^2.  Flipping s_i changes Q by
   sum_{e cut by i} (4 - 4 s_i delta_e(i) <s,delta_e>) <= 0, so with x_e := <s,delta_e>/u in [-1,1]:
   s_i sum_e x_e delta_e(i) >= M_i/u for all i.  A vertex of the polytope {x in [-1,1]^E: s_i sum_e x_e delta_e(i) >= M_i/u}
   has at most u fractional coordinates; set them to 0 (loss <= u per direction).  QED
   General r: choose signs sigma_i on the groups of pi_i (a cut of the quotient multigraph with >= M_i/2 wires between
   + and - groups), delta_e(i) = (sigma_i(head) - sigma_i(tail))/2, and b_i >= |sum_e x_e delta_e(i)|:
       log2 N >= min_i M_i/(2u) - u.
   With (3.2):  THEOREM E.  On a network of conductance phi, every DF identity whose u directions are beta-balanced
   has  N >= 2^{phi' M/(2u) - u}  (r = 2: 2^{phi' M/u - u}).   So on expanders "2^{M/u}" cannot be beaten by more than
   the constant phi'/2 in the exponent at one level with u directions -- for ANY labelling and ANY identity.
   Checked: c1 (263 random small labellings: u*A''/min M_i >= 1.000), c3 (expander, all cases >= 1.5).

(3.7) The flattening method CANNOT give M/sqrt(u) for arbitrary labellings [P + N]  (new, negative).
   Divergence form of A'' [P]: for an oriented selection x let y_v = (#selected wires with row-index end at v) - (#with
   column-index end at v).  Then e1_i(g) - e2_i(g) = y(g):  A'' = max_x min_i (1/2) sum_{g in pi_i} |y(g)|; the bound
   depends on x only through its divergence, |y(S)| <= cut(S).
   CIRCULAR LABELLING: 2u classes K_0..K_{2u-1} on a cycle, direction i (0 <= i < u) = the "diameter cut"
   H_i = K_i u ... u K_{i+u-1} versus its complement.  Every direction is a balanced cut (beta = 1/2 for equal classes).
   Prop. C.  For every network and every oriented selection, min_i b_i <= (1/2) max_i (cut(K_i) + cut(K_{i+u}));
   with classes of equal volume 2M/(2u) this is <= M/u.  So A'' <= M/u.
   Proof.  F_i := y(H_i), b_i = |F_i| (y(V) = 0).  F_{i+1} - F_i = y(K_{i+u}) - y(K_i) and F_u = y(V \ H_0) = -F_0.  The
   sequence F_0..F_u changes sign, so two consecutive terms have opposite signs and the smaller is <= half their
   difference; |F_u| = |F_0| is a direction too.  |y(K)| <= cut(K) <= vol(K).  QED
   Numerics (c3_expander.py, LP-exact up to +-u): random 4-regular graph, n = 48, M = 96, vertices assigned at random:
       product labelling {0,1}^u:  A''/M = 0.406, 0.340, 0.285  (u = 2,3,4)    ~ 0.564/sqrt(u) = 0.399, 0.326, 0.282
       circular labelling:         A''/M = 0.281, 0.240, 0.1667, 0.1250, 0.0833 (u = 3,4,6,8,12) = 1/u for u >= 6,
       while min_i M_i/M stays in [0.40, 0.50].   Same picture on a random 3-regular graph.
   Class-level LP (c2_worst.py, staircase = circular folded by the antipode): V = 1/u exactly for u = 2..8.
   Together with Lemma Q: min over labellings of A'' is Theta(min_i M_i/u), attained by the circular labelling.
   CONSEQUENCE.  A theorem "on expanders every labelling costs >= 2^{c M/sqrt u}" is NOT provable by A'' or any
   flattening (rank of one restriction): for the circular labelling every flattening certifies at most M/u.
   Moreover the circular labelling has no packing constraint (u may be as large as n/2), so by flattening alone one
   cannot even exclude N = 2^{O(sqrt M)} there.  This is a limit of the PROOF METHOD, not a construction:
   the only known identity on the circular labelling is cut-conditioning, 2^{M_i} >= 2^{phi' M} terms.
   Evidence that the truth is above M/u there [E]: the smallest staircase instance (u = 3, three antipodal wires
   (0,7),(1,6),(2,5), M = 3) has A'' = 1 but minimum N = 4 in the GF(2) sub-model (G18s table: N <= 3 UNSAT for all
   profiles) and N = 2, 3 not found numerically over R (G18e).  For u = 2 the staircase packet is the Hurwitz / HWY
   instance of G18e sec. 3 (N between n and O(n^2/log n), n = 2^{M/2}; a bound n^{1+eps} is a known open problem), so
   improving on flattening for such labellings is at least as hard as that problem [E, by analogy for u > 2].

## 4. Intermediate regime: A'' as a function of the load distribution over Hamming distances
EXCHANGEABLE loads (r = 2, or halves of [r]): a fraction p_d of the wires at Hamming distance d (d = 0..u), spread
uniformly over the sign types of that distance (which coordinates differ, with which orientations).  Then [P]
      A''/M = (1/u) sum_d p_d E|S_d|,    E|S_d| = d C(d-1, floor((d-1)/2))/2^{d-1} ~ sqrt(2d/pi)   (S_d = +-1 walk).
 Proof.  >=: select every wire with <1,delta_e> != 0, row-index end = the side with more +; a wire of distance d
 containing i contributes E|S_d|/d on average, and P(i in support) = d/u.  <=: for every x, s:
 min_i <= mean_i = (1/u) sum_e x_e <s,delta_e> <= (1/u) sum_e |<s,delta_e>|, and by exchangeability this is the same for
 all s.  QED     Checked by LP for u = 2..5, 30 load vectors (c4_symmetric.py: 30/30 equal to 1e-7).
 Values: d = 1 only: M/u (the Hamming-1 construction; equals min M_i).  d = u only: 0.80 M/sqrt(u).
 Product-uniform (d ~ Bin(u,1/2)): M C(2u,u)/4^u = 0.564 M/sqrt(u)  (G18s).
 In terms of the cut fraction m = M_i/M = mean(d)/u (same for all i by exchangeability):
      0.707 * m * M/sqrt(u)  <=  A''  <=  sqrt(m) * M/sqrt(u)   (E|S_d|/d decreasing, sqrt(d/2) <= E|S_d| <= sqrt(d); 0.80 asymptotically on both sides)
 On an expander m >= phi':  A'' >= 0.707 phi' M/sqrt(u).  The M/u value needs mean(d) = O(1), i.e. m = O(1/u), i.e. a
 direction that is a separator of size O(M/u): excluded as soon as u > 1/phi'.
 ANSWER for exchangeable loads: NO regime strictly between M/sqrt(u) and M/u on expanders; A'' = Theta_phi(M/sqrt u).
 Without exchangeability the formula is not a lower bound (circular labelling: formula 0.148, true A'' 0.083 at u = 12);
 there the proved range on an expander is  phi' M/(2u) - u <= log2 N <= min_i M_i, and A'' itself can be anywhere in
 [min M_i/u, min M_i] (both ends attained: circular, Hamming-1).
 For the A'' BOUND the whole interval between M/u and M/sqrt(u) is "reachable on hard instances" (circular labelling
 with part of the classes merged); for actual IDENTITIES nothing below 2^{phi' M} is known outside App. D's
 near-uniform packets.

## 5. Verdict: PARTLY
 CLOSED [P]: the concrete opening (loads on Hamming-distance-1 classes => N = 2^{min M_i} <= 2^{M/u}).  It is the
   cut-conditioning identity; on a network of conductance phi every balanced direction cuts >= phi' M wires, so
   (i) Hamming-1 concentration needs u <= 1/phi', (ii) more generally at most a fraction 1 - (phi' u - 1)/(u - 1) of the
   wires can be at distance <= 1, (iii) the identity has >= 2^{phi' M} terms, (iv) for exchangeable loads the lower
   bound stays >= 0.707 phi' M/sqrt(u), (v) for ARBITRARY balanced labellings and arbitrary DF identities
   N >= 2^{phi' M/(2u) - u} (Theorem E).  The designer's freedom in the labelling is real (Lemma D.3 is a design
   choice) but buys nothing through this identity on hard instances.
 NOT CLOSED: "every labelling of an expander costs >= 2^{c M/sqrt u}" -- unproved, and unprovable by flattening
   (circular labelling: A'' <= M/u on every network, Prop. C).  Between M/(2u) and M/sqrt(u) nothing is proved for
   non-exchangeable labellings; with u unconstrained (circular labelling needs no packets) Theorem E alone does not even
   give 2^{Omega(M/log M)}.
 NOT OPEN either: no identity with 2^{o(M)} terms is known on any non-uniform labelling of an expander.
 Precise open problem (smallest form): N_DF for the staircase packet -- u bundles of m wires, bundle c joining class c
   to class c+u on the 2u-cycle, directions = diameter cuts.  Flattening: N >= 2^m.  Trivial: N <= 2^{um}.
   u = 2 is the Hurwitz/HWY problem; u = 3, m = 1: N = 4 in the GF(2) sub-model.  A 2^{O(m)} identity for growing u
   would reopen the door (exponent M/u with u not limited by log M); a 2^{Omega(um/sqrt u)} bound needs a non-flattening
   technique.
 Unproved steps / hypotheses: (a) Theorem E is one-level, counts terms, needs beta-balanced directions and
   conductance (not treewidth) of the network actually contracted; (b) conductance of the 3-SAT group graph after
   App. D's encoding not computed; (c) App. D's own upper bound still taken at face value; (d) "truth above M/u on the
   circular labelling" rests on one 3-wire instance in a sub-model; (e) r >= 3: constant 1/2 lost (max-cut of quotient),
   and the cut-conditioning identity is not optimal (matrix-multiplication tensors appear).

## 6. Adversarial self-review
 * Lemma Q: the flip computation uses delta_e(i)^2 = 1 on cut wires and (a - 2b)^2 - a^2 = 4 - 4ab with b = s_i delta_e(i);
   x_e = <s,delta_e>/u is in [-1,1] because |<s,delta_e>| <= d_e <= u.  Rounding: basic solutions of a box LP with u
   extra constraints have <= u non-box-tight coordinates.  Checked numerically only on small/medium cases.
 * Prop. C uses only y(V) = 0 and |y(K)| <= cut(K); it does not need balanced classes for the first inequality.
   It bounds A'' (all flattenings "set unused bits to 0, one row end per wire"); restrictions that fix unused bits to
   other constants give the same matrices up to the same argument (I restricted is 0 unless both ends agree), but
   methods other than a single flattening rank (e.g. combining several flattenings non-linearly, Young flattenings,
   2-adic Smith profiles) are NOT covered.
 * (3.2) needs conductance in the volume sense; for bounded-degree graphs this is edge expansion / max degree.
 * The c2 local search did not find 1/u for u = 4, 5 (found 1/3, 0.26): the search is weak; the value 1/u is
   established by the explicit staircase (LP, u <= 8) and by Prop. C, not by the search.
 * c3 LP values are for one random graph per degree and one random assignment; A''_LP is within u of the integer A''.
 * Nothing here is a statement about 3-SAT or ETH; all bounds are in model DF (one level, term count).

Files (all in this folder): notes.md; c1_construct.py / c1_out.txt (exact identity + A'' brute force);
c2_worst.py / c2_pm_*.txt, c2_tern_3.txt (worst label measure search); c3_expander.py / c3_out.txt (expander LPs);
c4_symmetric.py / c4_out.txt (exchangeable-load formula); c5_staircase.py / c5_out.txt (staircase V = 1/u).
Note: the 3-regular run of c3 was stopped after circular u = 6 (time); c2 runs for u = 5, 6 were stopped early.
